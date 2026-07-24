"""
sprinx.label: Sprinzl-coordinate annotation logic for mt-tRNAs.

1. problem
   mt-tRNAs exist in four structural shapes: cloverleaf, D-armless, T-armless,
   doubly-armless (Ozerova et al. 2024). Sprinzl labels must be
   assigned relative to the correct shape; wrong shape -> wrong labels.

2. why score-based CM selection fails
   E-values are calibrated per model (Infernal User Guide);
   an armless CM with fewer columns produces better E-values for canonical
   sequences than the canonical CM does, regardless of biological fit.
   length-normalising (bits/column) doesn't help: armless CMs retain the
   highest-information columns (acceptor + anticodon stems), inflating
   per-column scores. Rfam avoids this with hand-set per-family GA cutoffs;
   this module avoids it by never comparing scores across models at all.

3. pipeline
   a. align to a canonical CM with cmalign --notrunc --nonbanded -g. --canonical-cm
      accepts multiple sources tried in priority order (e.g. bacterial whole-family
      CM, then a metazoan per-AA directory); the first tier whose anticodon anchors
      unambiguously and threads its own stem to the full canonical 5bp wins.
      a short thread never disqualifies a tier outright, since a real anticodon
      stem can be shorter than 5bp, so the best-threaded anchored tier
      is kept and used if no later tier does better. never by score/E-value
      (see (2)): a CM built for the wrong clade can fail to thread a divergent
      sequence at all, or mis-thread the otherwise-invariant anticodon stem
      specifically even when the anticodon itself anchors. details in
      select_cm_and_align.
   b. anchor on the anticodon; a missing UPSTREAM arm (D-arm) shifts remaining
      structure into wrong model columns (register shift). missing DOWNSTREAM
      arm (T-arm) does not shift. measure offset = expected_anticodon_slot - observed.
   c. n_pairs==0 at a stem slot means zero alignment columns have BOTH pairing
      partners simultaneously non-gap. no column can form a pair, so no stem can
      exist there: geometry forces the call, with no threshold to tune.
      n_pairs==0 has two distinct causes that require different responses:
        (i)  arm loss: the sequence simply has no arm. the element span
             across the alignment is mostly or entirely gap characters.
        (ii) CM threading failure: the arm exists but cmalign placed its sequence
             into unmodeled insert columns because the arm is too divergent from
             the CM consensus. the stem model columns are all gaps, but the span
             DOES contain nucleotides as insert characters.
      distinguishing (i) from (ii): count non-gap nucleotides across the full
      element span (stem + loop model columns + intervening insert characters).
      if the count is < n_stem_cols + MIN_HAIRPIN_LOOP (=3, steric minimum for
      the RNA backbone to close a hairpin), no hairpin can form physically:
      arm loss (i). otherwise: threading failure (ii).
      hybrid Infernal + RNAfold design: for threading failures, Infernal's
      canonical CM is correct for all arms it DID thread properly; only the
      mis-threaded arm needs structural recovery. RNAfold MFE on the short arm
      span (typically 13-20 nt) is reliable at this length because competing folds
      are energetically negligible. the hybrid avoids two failure modes: (a)
      relying on Infernal alone would call threading failures as arm loss and
      misroute to an armless CM; (b) relying on RNAfold alone for full-sequence
      mt-tRNA folding is unreliable due to tertiary interactions and base
      modifications not captured by 2D MFE.
   d. if truly absent, reroute to armless CM (Ozerova et al. 2024).
      isoacceptors (Leu1/Leu2, Ser1/Ser2) disambiguated by
      anticodon, not filename suffix. for doubly-armless (D + T both missing),
      routes to the d_and_t CM.
   e. assign Sprinzl coordinates.

4. implementation notes
   cmalign flags (required together, every call):
     --notrunc   : include all positions; without it, local mode silently drops
                   regions that fit poorly, causing false arm-loss calls.
     --nonbanded : exact CYK/Inside DP; HMM banding is ~10x faster but
                   introduces alignment errors on divergent mt-tRNA structures.
     -g          : glocal; prevents local begin/end states skipping arm regions.
   header format (pipe-delimited):
     field 1: seq id | field 2: three-letter aa (e.g. Ala, Leu1)
     field 3: anticodon (3nt, RNA or DNA) | field 4: taxon
     fallback 1: 'anticodon=XXX' tag anywhere in the header.
     fallback 2: GtRNAdb-style 'tRNA-{AA}-{anticodon}' name anywhere in the
     header (e.g. mt-tRNA-Ala-TGC-1-1); aa has no isoacceptor digit in this
     convention (Leu/Ser cover both isoacceptors), so aa_field_to_cm_code
     returns the bare code and CM selection disambiguates by anticodon anchor
     (_pick_by_anticodon_anchor), same as filename-suffixed isoacceptor CMs.
     field 3 (or the fallback anticodon) is the primary key for CM selection;
     field 2 (or the fallback aa) only identifies aa.
   armless CM filenames: armless_trn{AA}_wo_{arm}.cm where arm is d, t, or d_and_t
   for doubly-armless (Ozerova et al. 2024). armless CM rerouting
   is unaffected by which canonical CM tier won above; it only triggers once
   an arm-loss diagnosis is made from whichever tier's alignment was used.
   each --canonical-cm source is a directory of {label}_{AA}.cm files (e.g.
   Metazoan_P.cm; label/clade is ignored, selection is by AA only, per-sequence,
   same as armless CM selection) or a single CM file (applies to every aa,
   e.g. a whole-family CM like TRNAinf-bact.cm).

5. output
   sprinzl_mapping.tsv: seq_id, seq_index, nucleotide, sprinzl_position, region,
   cm_used, rerouted, arm_loss_call, structure (dot-bracket symbol at this
   position). the structure column lets scripts/visualize_ss.py reconstruct
   each record's secondary structure from the TSV alone, with no need to
   re-run cmalign; see scripts/visualize_ss.py for optional R2DT-rendered
   2D diagrams (needs its own extra dependencies and a Singularity/R2DT image,
   not required for the core sprinx package).
"""

import os
import re
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict

import RNA
from forgi.graph.bulge_graph import BulgeGraph
from Bio.Data.IUPACData import protein_letters_3to1
from loguru import logger

warnings.filterwarnings("ignore")

def _configure_logging(level):
    """(re)point loguru at stderr with a bare message format. called at import
    (INFO) and again wherever --debug is honoured (main(), each worker process
    since multiprocessing forks/spawns fresh interpreters)."""
    logger.remove()
    logger.add(sys.stderr, format="<level>{message}</level>", level=level)


_configure_logging("INFO")


# --- constants: tRNA topology facts + Sprinzl coordinate system ---

WC_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}

# anticodon arm is the 2nd inner stem-loop (0-indexed) in a canonical cloverleaf;
# topological fact, not tunable; changing it requires a different CM.
EXPECTED_ANTICODON_STEM_INDEX = 1

SPRINZL_REGION = {}
for _p in range(1, 8):   SPRINZL_REGION[str(_p)] = "acceptor_5"
for _p in (8, 9):         SPRINZL_REGION[str(_p)] = "connector_AD"
for _p in range(10, 14):  SPRINZL_REGION[str(_p)] = "D_stem_5"
for _p in list(range(14, 22)) + ["17a", "20a", "20b"]:
    SPRINZL_REGION[str(_p)] = "D_loop"
for _p in range(22, 26):  SPRINZL_REGION[str(_p)] = "D_stem_3"
SPRINZL_REGION["26"] = "connector_DC"
for _p in range(27, 32):  SPRINZL_REGION[str(_p)] = "C_stem_5"
for _p in range(32, 39):  SPRINZL_REGION[str(_p)] = "C_loop"
for _p in range(39, 44):  SPRINZL_REGION[str(_p)] = "C_stem_3"
for _p in (44, 45):       SPRINZL_REGION[str(_p)] = "V_loop"
for _p in range(46, 49):  SPRINZL_REGION[str(_p)] = "V_loop"
for _p in range(1, 8):    SPRINZL_REGION[f"e1{_p}"] = "V_stem_5"
for _p in range(1, 6):    SPRINZL_REGION[f"e{_p}"] = "V_loop"
for _p in range(1, 8):    SPRINZL_REGION[f"e2{_p}"] = "V_stem_3"
for _p in range(49, 54):  SPRINZL_REGION[str(_p)] = "T_stem_5"
for _p in range(54, 61):  SPRINZL_REGION[str(_p)] = "T_loop"
for _p in range(61, 66):  SPRINZL_REGION[str(_p)] = "T_stem_3"
for _p in range(66, 73):  SPRINZL_REGION[str(_p)] = "acceptor_3"
for _p in range(73, 77):  SPRINZL_REGION[str(_p)] = "discriminator_CCA"

# armless CM filename regex; naming follows Ozerova et al. 2024.
ARMLESS_CM_RE = re.compile(r"armless_trn(\w+)_wo_(d_and_t|d|t)\.cm$")


# --- generic helpers ---

def run(cmd):
    """run a subprocess; log the exact command at debug level for manual reproduction."""
    logger.debug(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.stdout, result.stderr, result.returncode


def drop_orphan_brackets(ss):
    """after stripping CM gap columns, a deletion on a paired column leaves its
    partner bracket dangling. converts unpartnered '(' or ')' to '.' so
    RNA.ptable() will accept the result. hand-rolled stack walk (not a library
    call) because this is specific to the gap-stripping context."""
    ss = list(ss)
    stack = []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                stack.pop()
            else:
                ss[i] = "."
    for i in stack:
        ss[i] = "."
    return "".join(ss)


HEADER_TRNA_NAME_RE = re.compile(r"tRNA-([A-Za-z]{3})-([ACGTUacgtu]{3})")


def header_to_anticodon(header):
    """extract anticodon from 'id|aa|anticodon|taxon' (field 3), 'anticodon=XXX'
    tag, or a GtRNAdb-style 'tRNA-{AA}-{anticodon}' name (e.g.
    'mt-tRNA-Ala-TGC-1-1') anywhere in the header. returns 3-nt RNA string or
    None; returns None rather than guessing on format mismatch, since a wrong
    anticodon propagates through the entire Sprinzl assignment."""
    fields = header.split("|")
    if len(fields) >= 3 and re.fullmatch(r"[ACGUTacgut]{3}", fields[2]):
        return fields[2].upper().replace("T", "U")
    m = re.search(r"anticodon=([ACGUTacgut]{3})", header)
    if m:
        return m.group(1).upper().replace("T", "U")
    m = HEADER_TRNA_NAME_RE.search(header)
    return m.group(2).upper().replace("T", "U") if m else None


def aa_field_to_cm_code(aa_field, cm_index_keys):
    """header aa field (e.g. 'Ala', 'Leu1') -> one-letter CM code (e.g. 'A', 'L1').
    derivation: strip digit suffix, protein_letters_3to1 (IUPAC 1984), reattach
    suffix, check against cm_index_keys. a bare aa field with no isoacceptor
    digit (e.g. 'Leu' from a GtRNAdb-style header, which never numbers
    isoacceptors) is returned as-is if it matches the digit-stripped form of
    one or more index entries, so the caller can disambiguate by anticodon
    (resolve_armless_cm, _resolve_canonical_for_tier) instead of failing here.
    returns None if the code matches nothing in the index at all."""
    if not aa_field:
        return None
    m = re.fullmatch(r"([A-Za-z]+)(\d*)", aa_field.strip())
    if not m:
        return None
    one = protein_letters_3to1.get(m.group(1).capitalize())
    if one is None:
        return None
    code = one + m.group(2)
    keys = {aa for aa, _ in cm_index_keys}
    if code in keys or any(k.rstrip("0123456789") == code for k in keys):
        return code
    return None


def header_to_aa(header):
    """return the aa field: the second pipe-delimited field, or (fallback) the
    aa name from a GtRNAdb-style 'tRNA-{AA}-{anticodon}' header. None if
    neither is present."""
    fields = header.split("|")
    if len(fields) >= 2 and fields[1].strip():
        return fields[1].strip()
    m = HEADER_TRNA_NAME_RE.search(header)
    return m.group(1) if m else None


def header_to_taxon(header):
    """return the raw taxon/species field (fourth pipe-delimited field) or None."""
    fields = header.split("|")
    return fields[3].strip() if len(fields) >= 4 and fields[3].strip() else None


# --- CM library ---

def find_cm_files(cm_dir):
    """recursively list all .cm files under cm_dir."""
    return [
        os.path.join(root, f)
        for root, _, files in os.walk(cm_dir)
        for f in files if f.endswith(".cm")
    ]


def _scan_cm_files(cm_dir, pattern, key_fn, kind, exclude=None, warn_on_conflict=False):
    """shared walk-and-regex-match skeleton for the two CM index builders below.
    key_fn(match) turns a regex match into the index key; files that don't match
    (or match `exclude`, used to keep armless CMs out of the canonical index)
    are skipped with a debug log rather than guessed at, since mis-binning here would
    silently route to the wrong model. warn_on_conflict logs when a later file
    overwrites an earlier one under the same key, since that silently prefers
    one file over another rather than erroring."""
    index = {}
    for path in find_cm_files(cm_dir):
        base = os.path.basename(path)
        if exclude and exclude.search(base):
            continue
        m = pattern.search(base)
        if not m:
            logger.debug(f"  not a {kind} CM by naming convention, skipping: {path}")
            continue
        key = key_fn(m)
        if warn_on_conflict and key in index and index[key] != path:
            logger.warning(f"multiple {kind} CMs map to {key!r}: "
                           f"using {path} (overriding {index[key]})")
        index[key] = path
    return index


def index_armless_cms(cm_dir):
    """scan cm_dir for armless_trn{AA}_wo_{arm}.cm files; return {(aa, arm): path}."""
    index = _scan_cm_files(cm_dir, ARMLESS_CM_RE, lambda m: (m.group(1), m.group(2)), "armless")
    logger.info(f"indexed {len(index)} armless CMs: {sorted(f'{aa}/{arm}' for aa, arm in index)}")
    return index


# canonical CM filename regex: {label}_{AA}.cm, e.g. Metazoan_P.cm; label
# (clade or any prefix) is ignored, only the AA code after the last "_" is used.
CANONICAL_CM_RE = re.compile(r"^.+_(\w+)\.cm$")


def index_canonical_cms(cm_dir):
    """scan cm_dir for {label}_{AA}.cm canonical CM files; return {aa_code: path}.
    label (e.g. clade) is ignored; selection is by AA only. armless CM files
    are excluded even though they also match the generic {x}_{y}.cm shape. if
    the same AA code appears under multiple files, the last one found wins and
    a warning is logged, since this silently picks one label/clade over another."""
    index = _scan_cm_files(cm_dir, CANONICAL_CM_RE, lambda m: m.group(1), "canonical",
                           exclude=ARMLESS_CM_RE, warn_on_conflict=True)
    logger.info(f"indexed {len(index)} canonical CMs by aa: {sorted(index)}")
    return index


def _resolve_canonical_for_tier(header, seq, tier):
    """resolve one canonical-CM tier to a concrete .cm path for this header, or
    None if the tier doesn't apply (a per-AA dict with no entry for this aa).
    a plain path string tier applies unconditionally (e.g. a whole-family CM
    like TRNAinf-bact.cm, which models every amino acid with one CM). a bare
    aa code that matches more than one same-base entry (e.g. a GtRNAdb-style
    header's 'Leu' resolving to both L1 and L2, since that naming never
    carries an isoacceptor digit) is disambiguated by anticodon anchor via
    _pick_by_anticodon_anchor, the same approach resolve_armless_cm uses."""
    if not isinstance(tier, dict):
        return tier
    aa_code = aa_field_to_cm_code(header_to_aa(header), {(aa, None) for aa in tier})
    if aa_code is None:
        return None
    if aa_code in tier:
        return tier[aa_code]
    candidates = [path for code, path in tier.items()
                  if code.rstrip("0123456789") == aa_code.rstrip("0123456789")]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return _pick_by_anticodon_anchor(header, seq, header_to_anticodon(header), candidates)


# --- cmalign: one call per (sequence, CM), gapped alignment kept for classify_arm_loss ---

def parse_multi_sto(path_or_text, from_text=False):
    """parse a (possibly multi-seq) Stockholm file into ({name: aligned_seq}, ss_cons).
    used by tests only (test_data_bundle.txt); production always calls cmalign_one."""
    seqs = defaultdict(str)
    ss = ""
    lines = path_or_text.splitlines() if from_text else open(path_or_text, encoding="utf-8")
    for line in lines:
        line = line.rstrip("\n")
        if not line or line.startswith(("//", "#=GS", "#=GR")):
            continue
        if line.startswith("#=GC SS_cons"):
            ss += line.split(None, 2)[-1]
        elif not line.startswith("#"):
            parts = line.split()
            if len(parts) >= 2:
                seqs[parts[0]] += parts[-1]
    return dict(seqs), ss


def cmalign_one(header, seq, cm_path):
    """align one sequence to one CM with cmalign --notrunc --nonbanded -g
    (see module docstring for flag rationale). returns dict with aligned_seq,
    ss_cons, raw_sto, cm_path; or None on failure. gapped alignment is retained
    because classify_arm_loss needs alignment-column coordinates."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(f">{header}\n{seq}\n")
        fa_path = fh.name
    try:
        stdout, stderr, rc = run(["cmalign", "--notrunc", "--nonbanded", "-g", cm_path, fa_path])
        if rc != 0:
            logger.warning(f"cmalign failed (rc={rc}) for {header} against {cm_path}:\n{stderr.strip()}")
            return None
        aligned_seq, ss_cons = "", ""
        query_name = header.split()[0]
        for line in stdout.splitlines():
            if line.startswith("#=GC SS_cons"):
                ss_cons += line.split()[-1]
            elif not line.startswith(("#", "//")) and line.strip():
                parts = line.split()
                if parts and parts[0] == query_name:
                    aligned_seq += parts[-1]
        if not aligned_seq or not ss_cons:
            aln_st = "empty" if not aligned_seq else "ok"
            ss_st  = "empty" if not ss_cons else "ok"
            logger.warning(
                f"cmalign produced no usable alignment for {header} "
                f"against {cm_path} (aligned_seq={aln_st}, ss_cons={ss_st})"
            )
            return None
        return {"aligned_seq": aligned_seq, "ss_cons": ss_cons,
                "raw_sto": stdout, "cm_path": cm_path}
    finally:
        os.unlink(fa_path)


def finalize_structure(alignment):
    """gapped cmalign alignment -> (ungapped_seq, ungapped_dotbracket) for Sprinzl
    numbering. strips gap symbols ('-' deletion, '.' insert gap); converts WUSS
    to dot-bracket. ss_cons's pair can have one side gapped in THIS sequence
    (e.g. arm loss); nulling only that side and stripping would let
    naive re-matching silently re-pair the orphan with an unrelated stem
    (observed corrupting the acceptor stem). fix: read pairing from the full
    consensus db via RNA.ptable before stripping, and null BOTH sides of any
    pair where either column is gapped here, so stripping can't create an
    orphan at all."""
    aligned_seq, ss_cons = alignment["aligned_seq"], alignment["ss_cons"]
    db = list(RNA.db_from_WUSS(ss_cons))
    pt = RNA.ptable("".join(db))
    for i in range(len(db)):
        partner = pt[i + 1] - 1
        if partner >= 0 and (aligned_seq[i] in "-." or aligned_seq[partner] in "-."):
            db[i] = "."
    pairs = [(s, d) for s, d in zip(aligned_seq, db) if s not in "-."]
    seq = "".join(s for s, _ in pairs).upper().replace("T", "U")
    ss = drop_orphan_brackets("".join(d for _, d in pairs))
    return seq, ss


# --- structural arm-loss detection
#
# replaces E-value CM selection. the core problem with gap-fraction per stem:
# a missing UPSTREAM arm (D-arm) causes cmalign to slide the remaining sequence
# one slot forward in the model (register shift), making the anticodon arm appear
# to occupy D-arm model columns. gap fraction at the anticodon slot is then zero
# not because the anticodon arm is missing but because cmalign put it elsewhere.
# fix: anchor on the anticodon from the header (never inferred), measure the
# offset between the expected stem index and where it actually landed.
# downstream arm loss (T-arm) does NOT shift the register; detected independently
# via n_pairs==0 at the T-arm slot. doubly-armless (D + T both missing) produces
# offset==0 with n_pairs==0 at both D-arm and T-arm slots simultaneously.

def _forgi_stem_groups(ss):
    """All physical stems via forgi's BulgeGraph.

    Why merging is needed: forgi splits one stem into two pieces wherever an
    interior-loop bulge interrupts it. Those pieces are merged back together
    here.

    How merge-vs-keep-separate is decided:
    - an interior-loop ('i') edge connects two stems whether it's a bulge
      inside one real helix, or a junction between two separate arms.
    - forgi's graph shape looks the same either way, so we count "anchors"
      instead: each hairpin is one anchor, plus the acceptor stem counts as
      one more anchor if the group also holds a hairpin.
    - more than one anchor: a real junction. keep the stems separate.
    - one or zero anchors: a bulge inside a single helix. merge.
    - the junction case only happens in armless/doubly-armless alignments,
      where two arms end up directly adjacent. a full cloverleaf's arms
      always join through a multiloop ('m'), never 'i', so they never merge.
    - see TestForgiStemGroups for the four cases this covers.

    Returns: list of dicts, sorted by span start, each with:
    - stem5_cols, stem3_cols, stem_cols (both sides)
    - loop_cols (empty for the acceptor stem - it has no hairpin loop)
    - span

    Shared by get_stem_loop_elements (arm-loss diagnosis) and parse_topology
    (Sprinzl labeling), so both use the same physical stems."""
    db = RNA.db_from_WUSS(ss)
    bg = BulgeGraph.from_dotbracket(db)
    stem_elems = sorted([e for e in bg.defines if e.startswith("s")],
                        key=lambda e: bg.defines[e][0])
    acceptor_elem = stem_elems[0] if stem_elems else None  # lowest start col

    def interior_neighbors(elem):
        return {nb2 for nb in bg.edges[elem] if nb.startswith("i")
                for nb2 in bg.edges[nb] if nb2.startswith("s") and nb2 != elem}

    def hairpin_neighbors(elem):
        return {nb for nb in bg.edges[elem] if nb.startswith("h")}

    def make_group(members):
        stem5 = sorted(col for g in members for col in range(bg.defines[g][0] - 1, bg.defines[g][1]))
        stem3 = sorted(col for g in members for col in range(bg.defines[g][2] - 1, bg.defines[g][3]))
        hairpins = {h for g in members for h in hairpin_neighbors(g)}
        loop_cols = []
        if hairpins:
            a, b = bg.defines[next(iter(hairpins))][:2]
            loop_cols = list(range(a - 1, b))
        stem_cols = sorted(stem5 + stem3)
        return {"stem5_cols": stem5, "stem3_cols": stem3, "stem_cols": stem_cols,
                "loop_cols": loop_cols, "span": (min(stem_cols), max(stem_cols) + 1)}

    groups, used = [], set()
    for s in stem_elems:
        if s in used:
            continue
        group, frontier = {s}, [s]
        while frontier:
            cur = frontier.pop()
            for nb in interior_neighbors(cur):
                if nb not in group:
                    group.add(nb)
                    frontier.append(nb)

        hairpins = {h for g in group for h in hairpin_neighbors(g)}
        anchors = len(hairpins) + (1 if hairpins and acceptor_elem in group else 0)
        if anchors > 1:
            for g in group:  # real junction (acceptor<->arm or arm<->arm), not a bulge
                used.add(g)
                groups.append(make_group({g}))
        else:
            used |= group
            groups.append(make_group(group))

    groups.sort(key=lambda g: g["span"][0])
    return groups


def get_stem_loop_elements(ss):
    """ordered inner stem-loop dicts, excluding the outer acceptor stem (no
    hairpin loop of its own). each dict: {'stem_cols', 'loop_cols', 'span'}."""
    return [g for g in _forgi_stem_groups(ss) if g["loop_cols"]]


def _widen_arm_span(ss_cons, elements, idx):
    """widen elem['span'] from the CM's own column boundary out to the full
    gap between neighboring stems (or the acceptor). a threading failure can
    leave real arm sequence in columns the CM called flanking linker instead
    of its own; folding only elem['span'] then misses base pairs that belong
    to the same stem (pombe mt-Cys's D-arm recovers 3bp folded narrow, 5bp folded
    wide, since the extra 2bp were sitting in the linker). only ever call this
    on an ALREADY-confirmed threading failure (see select_cm_and_align);
    using it for detection itself lets real armless sequences fold a
    spurious hairpin out of unrelated linker sequence."""
    groups = _forgi_stem_groups(ss_cons)
    acceptor = next(g for g in groups if not g["loop_cols"])
    start = (elements[idx - 1]["stem3_cols"][-1] + 1 if idx > 0
             else max(acceptor["stem5_cols"]) + 1)
    end = (elements[idx + 1]["stem5_cols"][0] if idx + 1 < len(elements)
           else min(acceptor["stem3_cols"]))
    return start, end


def find_anticodon_stem_index(aligned_seq, stem_loop_elements, anticodon):
    """search for anticodon within each stem-loop's hairpin-loop columns only.
    both '-' and '.' must be stripped from loop sequences together; filtering
    only '-' can leave insert-column junk that produces a spurious extra match,
    breaking the "exactly one loop" assumption.
    returns (index, method) or (None, reason_string) on no/ambiguous match."""
    if anticodon is None:
        return None, "no_anticodon_in_header"
    candidates = [
        i for i, elem in enumerate(stem_loop_elements)
        if elem["loop_cols"] and anticodon in "".join(
            aligned_seq[c] for c in elem["loop_cols"] if aligned_seq[c] not in "-."
        )
    ]
    if len(candidates) == 1:
        return candidates[0], "unique_loop_match"
    if not candidates:
        return None, "no_loop_match"
    return None, f"ambiguous_{len(candidates)}_loop_matches"


def stem_complementarity(aligned_seq, ss, elem):
    """WC/wobble pairing check for one stem element. n_pairs: columns where both
    partners are simultaneously non-gap (0 pairs = structurally impossible for
    a stem to exist there, not a threshold call). n_compatible: of those, WC or
    G-U wobble pairs; callers read per_stem_complementarity directly rather
    than a binary verdict. raw WUSS in ss is handled transparently by
    db_from_WUSS."""
    db = RNA.db_from_WUSS(ss)
    pt = RNA.ptable(db)
    pairs = []
    for c in elem["stem_cols"]:
        partner = pt[c + 1]
        if partner > c + 1:  # count each pair once from 5' side
            a, b = aligned_seq[c], aligned_seq[partner - 1]
            if a not in "-." and b not in "-.":
                pairs.append((a, b))
    n = len(pairs)
    k = sum((a, b) in WC_PAIRS for a, b in pairs)
    return {"n_pairs": n, "n_compatible": k}


def classify_arm_loss(header, aligned_seq, ss_cons,
                      expected_anticodon_index=EXPECTED_ANTICODON_STEM_INDEX):
    """Top-level structural diagnosis for one cmalign'd sequence: which arm is
    missing.

    - D-arm, when a register shift occurs: measured via the shift itself.
    - D-arm (no shift) or T-arm: measured via per-slot absent() (see its
      docstring).

    Always returns full diagnostics for every stem, even on ambiguous input.

    See TestCanonical36, TestTArmless, TestDArmless, TestBothArmlessMature
    for end-to-end validation against real alignments."""
    anticodon = header_to_anticodon(header)
    elements = get_stem_loop_elements(ss_cons)
    n = len(elements)
    idx, method = find_anticodon_stem_index(aligned_seq, elements, anticodon)
    per_stem = [stem_complementarity(aligned_seq, ss_cons, e) for e in elements]

    result = {
        "anticodon": anticodon, "n_stem_loops": n,
        "anticodon_stem_index": idx, "anticodon_search_method": method,
        "register_offset": None, "per_stem_complementarity": per_stem,
        "call": "UNRESOLVED", "missing_arm": None,
    }

    def absent(i):
        """a stem counts as present only if both hold: (1) enough non-gap
        sequence occupies its columns at all (n_pairs >= MIN_STEM_PAIRS), and
        (2) enough of that sequence is actually WC/wobble-paired
        (n_compatible >= MIN_COMPATIBLE_PAIRS), not just coincidental residues
        sitting in aligned columns. absent() is the negation of that AND, so
        it's an OR of the two negated conditions: failing either one alone is
        enough to call the arm absent."""
        stem = per_stem[i]
        return stem["n_pairs"] < MIN_STEM_PAIRS or stem["n_compatible"] < MIN_COMPATIBLE_PAIRS

    if idx is None:
        # anticodon didn't anchor (ambiguous AT-rich anticodon): fallback scan.
        # no directional register shift available; can flag absence but not D vs T.
        flags = [i for i in range(n) if absent(i)]
        result["call"] = f"UNANCHORED_fallback_structurally_absent={flags}"
        # if both the first slot (D-arm equivalent) and last slot (T-arm) are
        # structurally absent, call doubly-armless even without an anticodon anchor
        if n >= 2 and 0 in flags and (n - 1) in flags:
            result["missing_arm"] = "d_and_t"
        elif flags:
            result["missing_arm"] = "ambiguous"
        return result

    offset = expected_anticodon_index - idx
    result["register_offset"] = offset

    if offset == 0:
        d_arm_idx = idx - 1
        t_arm_idx = n - 1
        d_absent = d_arm_idx >= 0 and absent(d_arm_idx)
        t_absent = t_arm_idx > idx and absent(t_arm_idx)
        # slots strictly between the anticodon and T-arm (e.g. an optional
        # variable-arm stem some CMs model) are reported but never decisive
        # on their own, since no armless CM exists to reroute a variable-arm loss to.
        other_missing = [i for i in range(idx + 1, t_arm_idx) if absent(i)]

        if d_absent and t_absent:
            # doubly-armless tRNAs show offset==0 because with both arms absent
            # the anticodon arm still lands in the expected model columns (no
            # single-arm register shift occurs).
            result["call"] = f"BOTH_ARMS_MISSING_slots={[d_arm_idx, t_arm_idx]}"
            result["missing_arm"] = "d_and_t"
        elif d_absent:
            # D-arm absent but NO register shift: cmalign left the D-arm's own
            # model columns gapped in place instead of sliding structure forward.
            # the shift isn't universal (seen with CMs modeling more than the
            # canonical D/C/T trio), so this direct per-slot check is needed too,
            # not just the offset>0 branch below.
            result["call"] = f"UPSTREAM_ARM_MISSING_slot={d_arm_idx}"
            result["missing_arm"] = "d"
        elif t_absent:
            missing_downstream = other_missing + [t_arm_idx]
            result["call"] = f"T_OR_VAR_ARM_MISSING_slots={missing_downstream}"
            result["missing_arm"] = "t"
        elif other_missing:
            result["call"] = f"T_OR_VAR_ARM_MISSING_slots={other_missing}"
        else:
            result["call"] = "CANONICAL_NO_ARM_LOSS"
    elif offset > 0:
        result["call"] = f"UPSTREAM_ARM_MISSING_offset={offset}"
        result["missing_arm"] = "d"
    else:
        # negative offset: not observed in any test fixture. surfaces explicitly
        # rather than silently coercing into another branch.
        result["call"] = f"UNEXPECTED_NEGATIVE_OFFSET={offset}"

    return result


# --- CM routing: canonical-first, structural diagnosis, RNAfold cross-check ---

# steric minimum: RNA backbone cannot close a hairpin with < 3 unpaired nts.
# triloops are the smallest observed RNA hairpins.
MIN_HAIRPIN_LOOP = 3

# soft threshold used by classify_arm_loss's absent(): 1-2 coincidental base
# pairs are too few to nucleate a stable helix, so n_pairs<3 is weak evidence
# a stem is real (3 is empirically the smallest count with no false positives
# on the canonical-36 test set). UNLIKE MIN_HAIRPIN_LOOP (geometric certainty),
# this is a judgment call, not a certainty; it only flags *candidates* for
# arm loss; every candidate still has to pass the hard arm_span_has_enough_sequence
# check before any reroute happens.
MIN_STEM_PAIRS = 3

# a single WC/wobble pair can't stack into a helix on its own.
MIN_COMPATIBLE_PAIRS = 2

# full canonical anticodon-stem length; a real stem can thread fewer pairs
# than this, so it's descriptive, not a hard minimum.
ANTICODON_STEM_PAIRS = 5


def arm_span_has_enough_sequence(aligned_seq, elem):
    """first-stage (fast, hard) filter after a stem slot is flagged absent: does
    the span contain enough nucleotides to physically form a hairpin
    (n_stem_cols + MIN_HAIRPIN_LOOP, the steric minimum)? False here means
    definite real loss. True is not proof of a real arm, just not ruled out by
    volume alone; see arm_is_threading_failure for the required 2nd check."""
    start, end = elem["span"]
    n_nts = sum(1 for c in aligned_seq[start:end] if c not in "-.")
    return n_nts >= len(elem["stem_cols"]) + MIN_HAIRPIN_LOOP


def _arm_full_span_subseq_and_fold(aligned_seq, final_seq, elem):
    """extract the FULL non-gap span (matched + insert columns together) and
    fold it with RNAfold MFE. full span, not insert-only: a real arm's
    sequence can land in the slot's own matched columns too (e.g. human
    mt-Val's T-arm under TRNAinf-bact.cm), which are safe to fold over since
    this slot was already flagged absent (n_pairs below MIN_STEM_PAIRS), so
    there are no real base pairs there to protect. shared by arm_is_threading_failure
    (detect) and patch_threading_failure_arm (source the patch).
    returns (ungapped_positions, arm_ss), or (None, None) if too little
    sequence to fold (< MIN_HAIRPIN_LOOP + 2 nt)."""
    span_start, span_end = elem["span"]
    ungapped_positions, ungapped_idx = [], 0
    for gapped_idx, c in enumerate(aligned_seq):
        if c not in "-.":
            if span_start <= gapped_idx < span_end:
                ungapped_positions.append(ungapped_idx)
            ungapped_idx += 1

    if len(ungapped_positions) < MIN_HAIRPIN_LOOP + 2:
        return None, None

    arm_subseq = "".join(final_seq[p] for p in ungapped_positions)
    arm_ss, _ = RNA.fold_compound(arm_subseq).mfe()
    return ungapped_positions, arm_ss


def arm_is_threading_failure(aligned_seq, final_seq, elem):
    """second-stage check, run only after arm_span_has_enough_sequence passes:
    does the span actually fold as a hairpin? needed because a CM with wide
    insert-state capacity can pass the raw-count check on unrelated leftover
    sequence alone. True: real, recoverable arm. False: real loss despite
    passing the count check."""
    _, arm_ss = _arm_full_span_subseq_and_fold(aligned_seq, final_seq, elem)
    return arm_ss is not None and "(" in arm_ss


def patch_threading_failure_arm(header, aligned_seq, final_seq, final_ss, elem):
    """Recover arm structure for a confirmed CM threading failure. Call only
    after arm_is_threading_failure.

    What it does: splices elem's own RNAfold fold into final_ss, rather than
    refolding the whole molecule. Full-sequence RNAfold on a mt-tRNA is
    unreliable (tertiary contacts, modified bases); a short isolated span
    (13-20nt) is fine.

    Why the span is cleared before the fold is written in: cmalign's
    structure inside this span is untrustworthy (that's why
    arm_is_threading_failure fired). A weak leftover bracket (below
    MIN_STEM_PAIRS but not blank) must not block the patch meant to replace
    it. Seen on real data: mt-Cys under TRNAinf-bact.cm had a single leftover
    pair blocking a 3bp D-stem fold that agreed with it and would have simply
    extended it.

    A bracket outside the span whose partner falls inside it would dangle
    once the span is cleared, so that partner is cleared too.

    Returns: patched final_ss, or the original final_ss if there's no fold
    or the result is unbalanced (safety net)."""
    ungapped_positions, arm_ss = _arm_full_span_subseq_and_fold(aligned_seq, final_seq, elem)

    if ungapped_positions is None or "(" not in arm_ss:
        return final_ss

    span = set(ungapped_positions)
    pt = RNA.ptable(final_ss)
    ss_list = list(final_ss)
    overridden, cleared_external = [], []
    for pos in ungapped_positions:
        if ss_list[pos] != ".":
            overridden.append((pos, ss_list[pos]))
        partner = pt[pos + 1] - 1
        if partner >= 0 and partner not in span:
            cleared_external.append(partner)
            ss_list[partner] = "."   # would otherwise dangle once span is cleared
        ss_list[pos] = "."

    if overridden or cleared_external:
        old_span = "".join(final_ss[p] for p in ungapped_positions)
        logger.warning(
            f"{header}: RNAfold patch overrode cmalign's own structure: "
            f"span positions {ungapped_positions[0]}-{ungapped_positions[-1]}: "
            f"old={old_span!r} -> new={arm_ss!r}; "
            f"{len(overridden)} bracket(s) inside the span replaced "
            f"({overridden}); "
            + (f"{len(cleared_external)} bracket(s) outside the span also cleared "
               f"to avoid dangling ({cleared_external})" if cleared_external
               else "none outside the span affected")
        )

    for i, ungapped_pos in enumerate(ungapped_positions):
        if i < len(arm_ss):
            ss_list[ungapped_pos] = arm_ss[i]

    patched = "".join(ss_list)
    if patched.count("(") != patched.count(")"):
        logger.warning(
            f"arm patch produced unbalanced structure "
            f"({patched.count('(')} opens, {patched.count(')')} closes); keeping original"
        )
        return final_ss
    return patched


def _pick_by_anticodon_anchor(header, seq, anticodon, candidates):
    """given >=2 CM paths that could all serve the same (ambiguous) aa code,
    align to each and keep whichever anchors the header's anticodon in its
    own anticodon loop. the anticodon is the discriminating fact, not
    filename or dict-key suffix.

    - shared by resolve_armless_cm (Leu1/Leu2, Ser1/Ser2 filenames) and
      _resolve_canonical_for_tier (same ambiguity, but from a bare
      GtRNAdb-style aa field with no isoacceptor digit at all).
    - falls back to the first candidate, logged, if none anchor or no
      anticodon is known.

    caveat: when two isoacceptor CMs model an identical anticodon-stem
    shape, this can pick either one, and the choice is arbitrary. checked
    against an Ascaris Leu2 sequence where both L1/L2 armless CMs
    anchored the anticodon identically - the "wrong" pick still produced
    byte-identical final structure and Sprinzl labels. a surprising cm_used
    value isn't itself proof of a labeling bug; confirm the actual output
    differs before treating it as one."""
    if anticodon:
        for path in candidates:
            aln = cmalign_one(header, seq, path)
            if aln is None:
                continue
            elements = get_stem_loop_elements(aln["ss_cons"])
            if find_anticodon_stem_index(aln["aligned_seq"], elements, anticodon)[0] is not None:
                logger.debug(f"{header}: isoacceptor {anticodon} -> {os.path.basename(path)}")
                return path
    logger.warning(f"{header}: isoacceptor disambiguation failed among "
                   f"{[os.path.basename(p) for p in candidates]}; using first")
    return candidates[0]


def resolve_armless_cm(header, seq, aa_code, missing_arm, anticodon, armless_cm_index):
    """pick the correct armless CM for this sequence and arm type. if multiple
    isoacceptor CMs share the same base aa code (Leu1/Leu2, Ser1/Ser2), disambiguate
    by anticodon anchor (_pick_by_anticodon_anchor). returns None if no CM exists."""
    if aa_code is None:
        return None
    candidates = [
        path for (code, arm), path in armless_cm_index.items()
        if arm == missing_arm and code.rstrip("0123456789") == aa_code.rstrip("0123456789")
    ]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return _pick_by_anticodon_anchor(header, seq, anticodon, candidates)


def _routing_result(final_alignment, cm_used, diagnosis, rerouted=False, threading_failure_elem=None):
    """assemble the dict select_cm_and_align returns at each of its exit points,
    so the shape is defined once instead of copy-pasted per branch."""
    return {"final_alignment": final_alignment, "cm_used": cm_used, "diagnosis": diagnosis,
            "rerouted": rerouted, "threading_failure_elem": threading_failure_elem}


def select_cm_and_align(header, seq, canonical_cm_tiers, armless_cm_index):
    """Top-level CM selection for one sequence.

    1. Align against every canonical CM tier (e.g. bacterial whole-family CM,
       then a metazoan per-AA directory). Never by raw alignment score
       (module docstring section 2).
       - Among tiers that anchor the anticodon, pick the one with the
         highest total base-pairing evidence summed across all stems
         (per_stem_complementarity's n_pairs). Ties keep the earlier tier.
       - Example: mt-Cys, TRNAinf-bact.cm threads 3/5 anticodon pairs,
         Metazoa_C.cm threads 5/5 for the identical sequence - the fuller
         thread wins.
       - Skip tiers that don't apply to this aa, whose cmalign fails, or
         that never anchor the anticodon.
       - If no tier ever anchors, fall back to the first tier that aligned
         at all.
    2. No arm missing (or only the variable arm): return the canonical
       alignment as-is.
    3. D-arm (no register shift) or T-arm flagged absent:
       - Cross-check with arm_span_has_enough_sequence, then
         arm_is_threading_failure, before trusting the flag.
       - Passing both: patch via RNAfold instead of rerouting.
       - D-arm via register shift (offset>0) skips this cross-check and is
         trusted directly - its span holds sequence displaced by the shift
         itself, not the D-arm's own (or absent) content, so both checks
         would false-positive there.
    4. Real loss: reroute via resolve_armless_cm (anticodon-disambiguated).
       No matching armless CM: warn and keep the canonical alignment.

    canonical_cm_tiers: list of tiers (path, or {aa_code: path} dict), or a
    bare path/dict wrapped as a single tier.

    Returns: dict with final_alignment, cm_used, diagnosis, rerouted,
    threading_failure_elem."""
    if isinstance(canonical_cm_tiers, (str, dict)):
        canonical_cm_tiers = [canonical_cm_tiers]

    canonical_alignment = canonical_cm = diagnosis = best_total_pairs = None
    first_alignment = first_cm = first_diag = None  # ultimate fallback: no tier ever anchors
    for tier in canonical_cm_tiers:
        path = _resolve_canonical_for_tier(header, seq, tier)
        if path is None:
            logger.info(f"{header}: skipping a canonical CM tier: no CM for this amino acid there")
            continue
        aln = cmalign_one(header, seq, path)
        if aln is None:
            logger.info(f"{header}: skipping a canonical CM tier: alignment failed")
            continue
        diag = classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
        if first_alignment is None:                    # ultimate fallback if nothing ever anchors
            first_alignment, first_cm, first_diag = aln, path, diag
        idx = diag["anticodon_stem_index"]
        if idx is None:
            logger.warning(
                f"{header}: anticodon did not anchor cleanly against {path}, "
                f"evaluating remaining canonical CM tiers\n"
                f"  aligned_seq={aln['aligned_seq']}\n"
                f"  ss_cons={aln['ss_cons']}"
            )
            continue

        total_pairs = sum(stem["n_pairs"] for stem in diag["per_stem_complementarity"])
        logger.debug(f"{header}: tier {path} anchors anticodon, total stem pairs={total_pairs}")
        if canonical_alignment is None or total_pairs > best_total_pairs:
            canonical_alignment, canonical_cm, diagnosis, best_total_pairs = aln, path, diag, total_pairs

    if canonical_alignment is None:
        canonical_alignment, canonical_cm, diagnosis = first_alignment, first_cm, first_diag
    if canonical_alignment is None:
        return _routing_result(None, None, None)

    missing_arm = diagnosis["missing_arm"]

    if missing_arm not in ("d", "t", "d_and_t"):
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    # step 3 (see docstring). D-arm via register shift is trusted directly: its
    # span contains non-D-arm sequence placed there by the shift itself, so
    # both checks would false-positive on truly D-armless sequences.
    # D-arm via no-shift (offset==0) doesn't have that problem and gets the
    # same cross-check as T-arm.
    if missing_arm in ("t", "d"):
        elem = idx = elements = None
        if missing_arm == "d" and diagnosis["register_offset"] != 0:
            pass  # register-shift D-arm: trusted directly, no cross-check
        else:
            elements = get_stem_loop_elements(canonical_alignment["ss_cons"])
            if missing_arm == "t":
                idx = len(elements) - 1 if elements else None
            else:
                idx = diagnosis["anticodon_stem_index"] - 1
                idx = idx if 0 <= idx < len(elements) else None
            elem = elements[idx] if idx is not None else None

        # detection stays on elem's own narrow span; only the confirmed
        # patch's fold gets widened (see _widen_arm_span).
        if elem and arm_span_has_enough_sequence(canonical_alignment["aligned_seq"], elem):
            final_seq, _ = finalize_structure(canonical_alignment)
            if arm_is_threading_failure(canonical_alignment["aligned_seq"], final_seq, elem):
                widened = _widen_arm_span(canonical_alignment["ss_cons"], elements, idx)
                wide_elem = dict(elem, span=widened)
                logger.info(
                    f"{header}: CM diagnosed {missing_arm}-arm missing against {canonical_cm} "
                    f"({diagnosis['call']}) but the span folds as a real hairpin "
                    f"(CM threading failure, not real arm loss); patching via RNAfold\n"
                    f"  aligned_seq={canonical_alignment['aligned_seq']}\n"
                    f"  ss_cons={canonical_alignment['ss_cons']}"
                )
                return _routing_result(canonical_alignment, canonical_cm, diagnosis,
                                        threading_failure_elem=wide_elem)

    # arm loss: reroute
    logger.info(f"{header}: {missing_arm}-arm missing against {canonical_cm} "
               f"({diagnosis['call']}), looking for an armless CM to reroute to\n"
               f"  aligned_seq={canonical_alignment['aligned_seq']}\n"
               f"  ss_cons={canonical_alignment['ss_cons']}")
    aa_code = aa_field_to_cm_code(header_to_aa(header), armless_cm_index.keys())
    anticodon = header_to_anticodon(header)
    armless_path = resolve_armless_cm(header, seq, aa_code, missing_arm, anticodon, armless_cm_index)

    if armless_path is None:
        logger.warning(f"{header}: {missing_arm}-arm missing ({diagnosis['call']}) "
                       f"but no armless CM for aa_code={aa_code!r}; using canonical")
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    armless_alignment = cmalign_one(header, seq, armless_path)
    if armless_alignment is None:
        logger.warning(f"{header}: armless CM realignment failed ({armless_path}); "
                       f"falling back to canonical despite {missing_arm}-arm loss")
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    return _routing_result(armless_alignment, armless_path, diagnosis, rerouted=True)


# --- topology + Sprinzl assignment ---

def parse_topology(ss):
    """Split the structure into the acceptor stem (the one stem with no
    hairpin loop of its own) and the inner stems.

    - Uses _forgi_stem_groups, which correctly keeps a bulge or a
      finalize_structure-nulled pair inside a stem (e.g. the acceptor's own
      5' or 3' half) as part of that same stem, rather than treating it as
      a phantom extra stem.
    - Does not label D/C/T yet - that needs the anticodon, see
      locate_anticodon_stem."""
    groups = _forgi_stem_groups(ss)
    acceptor = next((g for g in groups if not g["loop_cols"]), None)
    if acceptor is None:
        raise ValueError(f"no acceptor stem (no hairpin-less stem) in: {ss!r}")
    inner_stems = [g for g in groups if g["loop_cols"]]

    return {
        "acceptor_5": acceptor["stem5_cols"], "acceptor_3": acceptor["stem3_cols"],
        "inner_stems": inner_stems,
        "trailer": [p for p in range(acceptor["stem3_cols"][-1] + 1, len(ss)) if ss[p] == "."],
    }


# D always precedes C; any variable arm and T always follow C. So for a
# fixed total stem-loop count (other than 2, which is ambiguous either way
# depending on which arm is missing), the anticodon arm's position is fixed
# by topology alone: 3 stems (D, C, T) -> C is the middle one; 4 stems
# (D, C, variable arm, T) -> C is the second one.
EXPECTED_ANTICODON_ARM_INDEX = {3: 1, 4: 1}


def locate_anticodon_stem(topo, ss, seq, anticodon, missing_arm=None):
    """Identify the C-stem (anticodon arm), D-stem (sibling before C that
    does not enclose it), and T-stem (sibling after C).

    - inner_stems: a list of _forgi_stem_groups dicts (stem5_cols,
      stem3_cols, loop_cols already known from forgi).
    - 'does not enclose' check: a D-armless pseudostem opens before C and
      closes after C. It must be excluded as a D-arm candidate, since it IS
      the enclosing pseudostem, not a D-arm at all. See
      TestSprinzlAssignment::test_d_armless_replacement_loop_gets_d_arm_labels.

    How the C-stem is found: by position (EXPECTED_ANTICODON_ARM_INDEX). The
    anticodon sequence can coincidentally appear in more than one loop, so
    position is the deciding fact, not loop content.
    - The anticodon search still runs, but only as a sanity check: does the
      position-derived C-stem's own loop actually contain it? A mismatch
      means the alignment itself is broken, and raises rather than
      mislabeling silently.
    - Exactly 2 stems remaining is the one shape position alone can't
      resolve (C sits first if D is the missing arm, second if T is). That
      case needs missing_arm, already established via the register-offset
      diagnosis on the canonical alignment (classify_arm_loss).
    - Any other stem count, or 2 stems with missing_arm unknown, raises
      rather than guessing."""
    inner_stems = topo["inner_stems"]
    n_stems = len(inner_stems)

    def direct_loop(idx, group):
        # loop positions exclusive of any nested stem's span
        lp = set(group["loop_cols"])
        own_close = group["stem3_cols"][0]
        for j, other in enumerate(inner_stems):
            if j != idx and group["stem5_cols"][-1] < other["stem5_cols"][0] < own_close:
                lp -= set(range(other["stem5_cols"][0], other["stem3_cols"][-1] + 1))
        return sorted(lp)

    def unpaired(a, b):
        return [p for p in range(a, b) if ss[p] == "."]

    def loop_contains_anticodon(idx):
        return ac and ac in "".join(seq[p] for p in direct_loop(idx, inner_stems[idx]))

    ac = (anticodon or "").upper().replace("T", "U")
    by_position = sorted(range(n_stems), key=lambda i: inner_stems[i]["stem5_cols"][0])

    c_idx = None
    if n_stems == 2 and missing_arm in ("d", "t"):
        c_idx = by_position[-1] if missing_arm == "t" else by_position[0]
    elif n_stems in EXPECTED_ANTICODON_ARM_INDEX:
        c_idx = by_position[EXPECTED_ANTICODON_ARM_INDEX[n_stems]]
    elif n_stems == 1:
        c_idx = by_position[0]

    if c_idx is not None and ac and not loop_contains_anticodon(c_idx):
        raise ValueError(
            f"anticodon {ac!r} not found in the expected C-stem's own loop "
            f"(stem {c_idx} of {n_stems}); alignment looks broken, refusing "
            f"to guess which stem is the anticodon arm."
        )
    if c_idx is None and n_stems:
        raise ValueError(
            f"cannot determine the anticodon arm's position for a {n_stems}-stem "
            f"structure (anticodon={ac!r}, missing_arm={missing_arm!r})."
        )
    c_stem = inner_stems[c_idx] if c_idx is not None else None

    d_stem, t_stem, v_stem = None, None, None
    # outermost (last) column of the c-stem 3' strand; stem3_cols[0] is the
    # INNERMOST column (adjacent to the loop); using it here was a real bug:
    # var_loop's boundary must start after the whole c-stem ends, not after
    # its first (innermost) column, or var_loop's own assign_slots call
    # silently overwrites the c-stem-3 columns between [0] and [-1] with wrong
    # (var-loop) labels; confirmed on real data (e.g. mt-Glu): columns that
    # should be Sprinzl 40-43 (c_stem3) were coming out as 44-47 (var_loop).
    c_close = c_stem["stem3_cols"][-1] if c_stem else None
    if c_stem:
        before = [g for g in inner_stems
                  if g is not c_stem and g["stem5_cols"][0] < c_stem["stem5_cols"][0]
                  and g["stem3_cols"][0] < c_close]
        after = [g for g in inner_stems
                 if g is not c_stem and g["stem5_cols"][0] > c_close]
        d_stem = max(before, key=lambda g: len(g["stem_cols"])) if before else None
        # t-arm is the last (highest-position) stem after c-close.
        # min(after) breaks for class-ii tRNAs (ser, leu) and some mt-tRNAs
        # with a variable arm stem: it picks the variable arm as t-arm instead.
        t_stem = max(after, key=lambda g: g["stem5_cols"][0]) if after else None
        # a real variable-ARM stem (class-ii: Leu, Ser) is whatever's left
        # in `after` besides t_stem; only trusted when exactly one such
        # candidate remains; with a single "after" candidate there's no way
        # to tell a bare variable arm from a missing T-arm from topology
        # alone, so that ambiguous case is left to the existing missing_arm
        # machinery (classify_arm_loss) rather than guessed at here.
        v_candidates = [g for g in after if g is not t_stem]
        v_stem = v_candidates[0] if len(v_candidates) == 1 else None

    # outermost D-stem 3' column (strand edge), so the connector (pos 26) starts
    # only after the whole D-stem 3' strand; a D-stem-internal 3' bulge sits
    # before this edge and belongs to the D-stem, same strand-boundary reasoning
    # as the acceptor case below. using the innermost column instead would sweep
    # that bulge into linker_dc and label it 26 ahead of the real 24/25.
    d_stem3_end = d_stem["stem3_cols"][-1] if d_stem else None
    t_open = t_stem["stem5_cols"][0] if t_stem else None

    # acceptor STRAND boundaries, not member sets: an acceptor-internal bulge
    # (a '.' between two paired acceptor columns) belongs to the acceptor, so
    # var_loop/linker_5 must stop at the strand edge, not sieve out only the
    # literal paired columns and thereby swallow the bulge (which then gets a
    # wrong D-connector or V-loop label instead of an acceptor insertion).
    acceptor_3_start = topo["acceptor_3"][0]   # stem3_cols is sorted ascending
    acceptor_5_end = topo["acceptor_5"][-1]

    var_end = t_open if t_open is not None else (acceptor_3_start if c_close is not None else None)

    # class-ii variable ARM (Leu, Ser): a real nested stem-loop between the
    # c-arm and t-arm gets the Sprinzl e-series (e11-e17/e1-e5/e21-e27), not
    # the plain 44-48 sequential run; see sprinzl_map. ct_linker/vt_linker
    # are the (up to 2 / up to 3) unpaired nt flanking the v-stem
    # on either side; var_loop is the no-v-stem fallback (class-i short loop,
    # or no stem detected), using the previous flat 44-48 behaviour unchanged.
    v_stem5 = v_stem["stem5_cols"] if v_stem else []
    v_loop = v_stem["loop_cols"] if v_stem else []
    v_stem3 = v_stem["stem3_cols"] if v_stem else []
    if v_stem and c_close is not None and var_end is not None:
        ct_linker = list(range(c_close + 1, v_stem5[0]))
        vt_linker = list(range(v_stem3[-1] + 1, var_end))
        var_loop = []
    else:
        ct_linker = vt_linker = []
        if c_close is not None and var_end is not None:
            # all positions between c-stem close and t-arm open, not filtered to
            # unpaired: a variable arm stem's paired positions (when not split
            # out above) must still receive sprinzl 44-48 labels.
            var_loop = list(range(c_close + 1, var_end))
        else:
            var_loop = []

    d_stem5 = d_stem["stem5_cols"] if d_stem else []
    c_stem5 = c_stem["stem5_cols"] if c_stem else []
    t_stem5 = t_stem["stem5_cols"] if t_stem else []
    linker_5_end = d_stem5[0] if d_stem else (c_stem5[0] if c_stem else acceptor_3_start)

    # a short or RNAfold-patched T-stem can leave unpaired nts between its 3'
    # end and the acceptor 3' strand; fold them into the T-stem-3' run so
    # assign_slots numbers them (65, then 65A...) instead of leaving them blank.
    t_stem3 = t_stem["stem3_cols"] if t_stem else []
    if t_stem3:
        t_stem3 = t_stem3 + [p for p in range(t_stem3[-1] + 1, acceptor_3_start) if ss[p] == "."]

    return {
        "d_stem5": d_stem5, "d_stem3": d_stem["stem3_cols"] if d_stem else [],
        "d_loop": d_stem["loop_cols"] if d_stem else [],
        "c_stem5": c_stem5, "c_stem3": c_stem["stem3_cols"] if c_stem else [],
        "c_loop": c_stem["loop_cols"] if c_stem else [],
        "t_stem5": t_stem5, "t_stem3": t_stem3,
        "t_loop": t_stem["loop_cols"] if t_stem else [],
        "var_loop": var_loop,
        "ct_linker": ct_linker, "vt_linker": vt_linker,
        "v_stem5": v_stem5, "v_loop": v_loop, "v_stem3": v_stem3,
        # starts AFTER the acceptor's 5' strand ends: an acceptor-internal 5'
        # bulge sits before that edge and belongs to the acceptor, not the
        # linker (which otherwise mislabels it as a D-arm connector 8/9, or in
        # the D-armless case as a replacement-loop 8-26 position).
        "linker_5": unpaired(acceptor_5_end + 1, linker_5_end),
        "linker_dc": unpaired(d_stem3_end + 1, c_stem5[0]) if (d_stem3_end is not None and c_stem) else [],
    }


def assign_slots(labels, positions, slots):
    """map sequence positions onto Sprinzl slots 5'->3'. extras beyond len(slots)
    get letter-suffixed labels (e.g. '60A', '60B') on the last slot used."""
    for pos, slot in zip(positions, slots):
        labels[pos] = slot
    if len(positions) > len(slots) and slots:
        anchor = slots[-1]
        for k, pos in enumerate(positions[len(slots):]):
            letter = chr(ord("A") + k) if k < 26 else f"A{chr(ord('A') + k - 26)}"
            labels[pos] = f"{anchor}{letter}"


def _assign_anticodon_loop(labels, seq, c_loop, anticodon):
    """assign Sprinzl 32-38 anchored on the anticodon's own known position
    within the loop, not sequential order from the loop's 5' edge. a loop
    that isn't exactly the canonical 7nt (a stem-edge nucleotide the CM
    mis-threaded into the loop, or any other irregular split) would otherwise
    shift the anticodon off 34-35-36, the tool's core deliverable, even
    though its actual position is already known. among multiple coincidental
    matches within the loop, picks whichever is closest to the loop's own
    center, since the true anticodon is always the centered occurrence.
    fewer than 2 nt on either side of the anticodon leaves those canonical
    slots (32/33 or 37/38) simply unused, same as assign_slots elsewhere;
    more than 2 overflows via assign_slots' normal letter-suffix mechanism.
    falls back to sequential assignment only if the anticodon can't be
    located in the loop at all (e.g. no anticodon in header)."""
    ac = (anticodon or "").upper().replace("T", "U")
    loop_seq = "".join(seq[p] for p in c_loop)
    matches = [m.start() for m in re.finditer(f"(?={re.escape(ac)})", loop_seq)] if ac else []
    if not matches:
        assign_slots(labels, c_loop, [str(i) for i in range(32, 39)])
        return

    center = (len(loop_seq) - 3) / 2
    ac_start = min(matches, key=lambda i: abs(i - center))
    before, ac_positions, after = (
        c_loop[:ac_start], c_loop[ac_start:ac_start + 3], c_loop[ac_start + 3:])

    before_slots = ([str(i) for i in range(34 - len(before), 34)]
                     if len(before) <= 2 else ["32", "33"])
    assign_slots(labels, before, before_slots)
    for i, pos in enumerate(ac_positions):
        labels[pos] = str(34 + i)
    assign_slots(labels, after, ["37", "38"])


def sprinzl_map(ss, seq, anticodon, missing_arm=None):
    """Assign a Sprinzl label to every nucleotide index.

    - Returns: {seq_index: label}.
    - D-armless tRNAs: the replacement loop (all of linker_5) is mapped onto
      D-arm Sprinzl positions 8-26 by structural analogy.
    - Missing T-arm: no labels are produced for its region.
    - missing_arm (from classify_arm_loss's diagnosis on the canonical
      alignment, if this sequence was rerouted) is passed through to
      locate_anticodon_stem to resolve its 2-stem case; see its docstring."""
    topo = parse_topology(ss)
    arms = locate_anticodon_stem(topo, ss, seq, anticodon, missing_arm)
    labels = {}
    assign_slots(labels, topo["acceptor_5"], [str(i) for i in range(1, 8)])
    assign_slots(labels, topo["acceptor_3"], [str(i) for i in range(66, 73)])
    assign_slots(labels, topo["trailer"],    ["73", "74", "75", "76"])
    if arms["d_stem5"]:
        # canonical d-arm: stem/loop/connector assigned to their proper slots
        assign_slots(labels, arms["linker_5"],  ["8", "9"])
        assign_slots(labels, arms["d_stem5"],   [str(i) for i in range(10, 14)])
        assign_slots(labels, arms["d_loop"],
                     ["14", "15", "16", "17", "17a", "18", "19", "20", "20a", "20b", "21"])
        assign_slots(labels, arms["d_stem3"],   [str(i) for i in range(22, 26)])
        assign_slots(labels, arms["linker_dc"], ["26"])
    else:
        # d-armless: replacement loop occupies positions 8-26 sequentially.
        # locate_anticodon_stem puts all linker nucleotides into linker_5 when
        # d_stem is absent (linker_dc and d_stem3 are both empty in that case).
        _d_arm_slots = [
            "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "17a", "18", "19", "20", "20a", "20b", "21",
            "22", "23", "24", "25", "26",
        ]
        assign_slots(labels, arms["linker_5"], _d_arm_slots)
    assign_slots(labels, arms["c_stem5"],    [str(i) for i in range(27, 32)])
    _assign_anticodon_loop(labels, seq, arms["c_loop"], anticodon)
    assign_slots(labels, arms["c_stem3"],    [str(i) for i in range(39, 44)])
    if arms["v_stem5"]:
        # class-ii variable ARM (Leu, Ser): real nested stem-loop, Sprinzl
        # e-series; see locate_anticodon_stem's v_stem docstring. reserved
        # space is 7bp/5nt/7bp; anything beyond overflows via assign_slots'
        # usual letter-suffix mechanism (e.g. 'e17A'), same as every other
        # insertion-code slot in this scheme.
        # right-aligned on 45 (immediately before e11), same reasoning as
        # _assign_anticodon_loop's "before" segment: a single linker nt here
        # is adjacent to the stem and must get 45, not 44.
        ct_slots = ([str(i) for i in range(46 - len(arms["ct_linker"]), 46)]
                    if len(arms["ct_linker"]) <= 2 else ["44", "45"])
        assign_slots(labels, arms["ct_linker"], ct_slots)
        assign_slots(labels, arms["v_stem5"],   [f"e1{i}" for i in range(1, 8)])
        assign_slots(labels, arms["v_loop"],    [f"e{i}" for i in range(1, 6)])
        n3 = min(len(arms["v_stem3"]), 7)
        assign_slots(labels, arms["v_stem3"],   [f"e2{k}" for k in range(n3, 0, -1)])
        assign_slots(labels, arms["vt_linker"], ["46", "47", "48"])
    elif len(arms["var_loop"]) > 5:
        # no paired v-stem was threaded, but the variable region is longer
        # than the plain 44-48 5-slot capacity can hold; that excess length
        # is itself the signal of a real extended variable region (class-ii),
        # whether or not a stem happens to be threaded there. treat it as
        # loop-only: 44/45 before, e1-e5 (+letter-suffix overflow on e5 for
        # anything beyond 5nt) in the middle, 46/47/48 after, not a bigger
        # 44-48-derived overflow run, which was never a real Sprinzl code.
        var_loop = arms["var_loop"]
        before, middle, after = var_loop[:2], var_loop[2:-3], var_loop[-3:]
        assign_slots(labels, before, ["44", "45"])
        assign_slots(labels, middle, [f"e{i}" for i in range(1, 6)])
        assign_slots(labels, after,  ["46", "47", "48"])
    else:
        assign_slots(labels, arms["var_loop"], ["44", "45", "46", "47", "48"])
    assign_slots(labels, arms["t_stem5"],    [str(i) for i in range(49, 54)])
    assign_slots(labels, arms["t_loop"],     [str(i) for i in range(54, 61)])
    assign_slots(labels, arms["t_stem3"],    [str(i) for i in range(61, 66)])

    # strand ranges (first..last paired column of each stem strand); a
    # single-sided bulge is an unpaired column WITHIN one of these, which forgi
    # leaves outside the stem's own stem/loop columns. pass them so
    # _fill_stem_bulges only fills positions a stem actually owns.
    strands = [topo["acceptor_5"], topo["acceptor_3"]]
    for g in topo["inner_stems"]:
        strands += [g["stem5_cols"], g["stem3_cols"]]
    _fill_stem_bulges(labels, ss, strands)
    return labels


def _fill_stem_bulges(labels, ss, strands):
    """a single-sided bulge inside a stem (see _forgi_stem_groups: the bulged
    nucleotide sits outside the merged stem's own stem/loop columns, by
    construction) is a real nucleotide with no named Sprinzl slot; letter-
    suffix it onto the preceding assigned position, the same insertion-code
    convention assign_slots uses for trailing overhangs (60A, ...). ownership
    is checked against `strands` (each stem strand's first..last column span):
    only an unlabeled '.' that a stem actually spans counts as a bulge. every
    such bulge in real mt-tRNA data is a cmalign insert column (#=GC RF == '.'),
    never a model-consensus position, and the Sprinzl scheme has no canonical
    number for a stem bulge, so a suffix is the right label. mutates labels in
    place. an unlabeled '.' outside every stem strand, or any unlabeled paired
    column, is left alone so a different bug surfaces here rather than getting
    silently patched over."""
    owned = set()
    for strand in strands:
        if strand:
            owned |= set(range(min(strand), max(strand) + 1))
    last_label, next_ord = None, {}
    for i, c in enumerate(ss):
        if i in labels:
            last_label = labels[i]
        elif c == "." and i in owned and last_label is not None:
            n = next_ord.get(last_label, 0)
            next_ord[last_label] = n + 1
            letter = chr(ord("A") + n) if n < 26 else f"A{chr(ord('A') + n - 26)}"
            labels[i] = f"{last_label}{letter}"


# --- per-sequence worker: bundled for multiprocessing.Pool.map ---

def process_one_record(args):
    """worker for one (header, seq) FASTA record.

    - takes a single tuple for Pool.map compatibility.
    - canonical_cm_tiers and armless_cm_index are inside that tuple: each
      worker is a fresh process, and module-level globals aren't reliably
      shared across fork vs spawn.
    - per-tier canonical CM resolution (which .cm path applies to this aa,
      if any) happens inside select_cm_and_align; see
      _resolve_canonical_for_tier."""
    header, seq, canonical_cm_tiers, armless_cm_index, debug = args
    seq = seq.upper().replace("T", "U")

    if debug:
        _configure_logging("DEBUG")

    routing = select_cm_and_align(header, seq, canonical_cm_tiers, armless_cm_index)
    alignment = routing["final_alignment"]

    if alignment is None:
        logger.warning(f"{header}: cmalign failed, skipped")
        return {"header": header, "rows": [], "summary": "CMALIGN_FAILED"}

    final_seq, final_ss = finalize_structure(alignment)
    cm_only_ss = rnafold_only_ss = None
    if routing.get("threading_failure_elem"):
        cm_only_ss = final_ss
        # naive whole-sequence MFE fold, for comparison only; never used for
        # the actual patch (see module docstring: unreliable at full mt-tRNA
        # length, tertiary contacts and modified bases aren't 2D-foldable).
        rnafold_only_ss, _ = RNA.fold_compound(final_seq).mfe()
        final_ss = patch_threading_failure_arm(
            header, alignment["aligned_seq"], final_seq, final_ss,
            routing["threading_failure_elem"]
        )

    if len(final_seq) != len(seq):
        logger.warning(f"{header}: ungapped length {len(final_seq)} != input {len(seq)}, skipped")
        return {"header": header, "rows": [], "summary": "LENGTH_MISMATCH"}

    anticodon = header_to_anticodon(header)
    if anticodon is None:
        logger.warning(f"{header}: no anticodon in header; C-stem location unreliable")

    diagnosis = routing["diagnosis"] or {}
    sprinzl = sprinzl_map(final_ss, final_seq, anticodon, diagnosis.get("missing_arm"))

    unlabeled = [i for i in range(len(final_seq)) if i not in sprinzl]
    if unlabeled:
        logger.warning(f"{header}: {len(unlabeled)} position(s) left without a Sprinzl "
                       f"number at seq index {unlabeled}; output rows for them are blank")

    cm_name = routing["cm_used"] or "NONE"
    if cm_name not in ("RNAfold", "NONE"):
        cm_name = os.path.basename(cm_name)
    summary = f"CM:{cm_name}" + (" [rerouted]" if routing["rerouted"] else "")

    logger.debug(f"{header}")
    logger.debug(f"  seq ({len(final_seq)}nt): {final_seq}")
    logger.debug(f"  ss  ({len(final_ss)}nt):  {final_ss}")
    logger.debug(f"  arm-loss: {diagnosis.get('call')}  "
                 f"(anchor:{diagnosis.get('anticodon_search_method')}, "
                 f"offset={diagnosis.get('register_offset')})")
    for i, stem in enumerate(diagnosis.get("per_stem_complementarity", [])):
        logger.debug(f"    stem[{i}]: n_pairs={stem['n_pairs']} "
                     f"n_compatible={stem['n_compatible']}")
    if routing.get("threading_failure_elem"):
        logger.debug(f"  threading failure: patched via RNAfold; ss: {final_ss}")
    if routing["rerouted"]:
        logger.debug(f"  rerouted to: {routing['cm_used']}")
    logger.debug(f"  raw stockholm:\n{alignment['raw_sto']}")

    rows = []
    for i, base in enumerate(final_seq):
        label = sprinzl.get(i, "")
        # 'e'-prefixed variable-arm labels (e11, e1, e23, ...) aren't purely
        # numeric like every other slot; match the optional 'e' along with
        # the digits so an overflow-suffixed one (e17A) still resolves to its
        # base code (e17) for the region lookup, same as '60A' -> '60' does.
        region_key = re.match(r"e?\d+", label).group() if label else ""
        rows.append({
            "seq_id": header, "seq_index": i, "nucleotide": base,
            "sprinzl_position": label, "region": SPRINZL_REGION.get(region_key, ""),
            "cm_used": cm_name, "rerouted": routing["rerouted"],
            "arm_loss_call": diagnosis.get("call"),
            # dot-bracket symbol at this position; carries final_ss into the
            # TSV so scripts/visualize_ss.py can rebuild structure per record
            # without re-running cmalign (see module docstring, section 5).
            "structure": final_ss[i],
            # pre-patch CM structure and naive whole-sequence RNAfold structure,
            # same indices as final_ss (patching replaces characters in place,
            # never changes length); blank when this record wasn't a threading
            # failure, matching cm_only_ss/rnafold_only_ss being None there.
            "cm_only_structure": cm_only_ss[i] if cm_only_ss else "",
            "rnafold_only_structure": rnafold_only_ss[i] if rnafold_only_ss else "",
        })

    logger.info(f"{header}: {summary}  [{diagnosis.get('call')}]")
    return {"header": header, "rows": rows, "summary": summary,
            "seq": final_seq, "ss": final_ss, "sprinzl": sprinzl,
            "cm_only_ss": cm_only_ss, "rnafold_only_ss": rnafold_only_ss}


# =============================================================================
# arm-loss call string glossary
# every processed sequence emits exactly one of these in the log (see
# README.md's "Arm-loss call glossary" for the same list, user-facing).
#
# CANONICAL_NO_ARM_LOSS: every stem-loop slot passes absent() (see classify_arm_loss).
#
# T_OR_VAR_ARM_MISSING_slots=[n,..]: one or more slots fail absent() (0-indexed,
#   5'->3'). a middle slot means an optional variable-arm stem, never decisive on
#   its own (no armless CM for it). the LAST slot is the T-arm: real loss, or a
#   threading failure patched via RNAfold; see select_cm_and_align step 3.
#
# UPSTREAM_ARM_MISSING_offset=n: D-arm absent, detected via register shift
#   (anticodon landed n slots downstream of expected). trusted directly, no
#   threading-failure cross-check (see select_cm_and_align).
#
# UPSTREAM_ARM_MISSING_slot=n: D-arm absent with NO register shift (its own
#   slot n=idx-1 fails absent() instead). seen with CMs modeling more than the
#   canonical D/C/T trio. DOES get the threading-failure cross-check, unlike the
#   offset-based call above.
#
# BOTH_ARMS_MISSING_slots=[n,..]: doubly-armless: D-arm and T-arm slots both
#   fail absent(), offset==0 (no single-arm shift with both arms gone). reroutes
#   to armless_trn{AA}_wo_d_and_t.cm.
#
# UNANCHORED_fallback_structurally_absent=[n,..]: anticodon not uniquely
#   anchored (ambiguous AT-rich triplet); no directional signal, less reliable.
#
# threading failure (separate log line, not a call string): "CM diagnosed X-arm
#   missing (...) but the span folds as a real hairpin ... patching via RNAfold".
#   patch aborts silently (DEBUG log) on a bracket conflict with existing structure.
# =============================================================================

