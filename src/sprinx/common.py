"""
sprinx.common: structural parsing and Sprinzl-label assignment shared by
sprinx.mito and sprinx.cyto.

This generic module (not mito or cyto-specific)
turns a cmalign alignment into stem/loop topology (via forgi) and maps that
topology onto Sprinzl coordinates, which both the mito and cyto paths
need identically.

cmalign flags (required together, every call):
  --notrunc   : include all positions; without it, local mode silently drops
                regions that fit poorly, causing false arm-loss calls.
  --nonbanded : exact CYK/Inside DP; HMM banding is ~10x faster but
                introduces alignment errors on divergent tRNA structures.
  -g          : glocal; prevents local begin/end states skipping arm regions.

header format (pipe-delimited):
  field 1: seq id | field 2: three-letter aa (e.g. Ala, Leu1)
  field 3: anticodon (3nt, RNA or DNA) | field 4: taxon
  fallback 1: 'anticodon=XXX' tag anywhere in the header.
  fallback 2: GtRNAdb-style 'tRNA-{AA}-{anticodon}' name anywhere in the
  header (e.g. mt-tRNA-Ala-TGC-1-1); aa has no isoacceptor digit in this
  convention (Leu/Ser cover both isoacceptors), so aa_field_to_cm_code
  returns the bare code and CM selection disambiguates by anticodon anchor
  (_pick_by_anticodon_anchor).
  field 3 (or the fallback anticodon) is the primary key for CM selection;
  field 2 (or the fallback aa) only identifies aa.
"""

import importlib.resources
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


HEADER_TRNA_NAME_RE = re.compile(r"tRNA-([A-Za-z]+\d*)-([ACGTUacgtu]{3})")


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
    instead of failing here. returns None if the code matches nothing in the
    index at all."""
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

def package_data_path(*parts):
    """resolve a path under sprinx's installed package data (src/sprinx/data/),
    e.g. package_data_path("cyto_cm", "TRNAinf-euk-iso"). used for the bundled
    default CM databases; assumes a normal (non-zipped) install."""
    return str(importlib.resources.files("sprinx").joinpath("data", *parts))


def check_cm_format(cm_path):
    """run cmstat on a CM file (single-model or multi-model) and raise a
    clear error if it fails. catches an unsupported format up front, at
    startup - e.g. old INFERNAL-1.0 CMs, which cmalign also refuses, but
    only after failing deep inside a worker process with a bare Infernal
    error and no indication which supplied CM caused it."""
    stdout, stderr, rc = run(["cmstat", cm_path])
    if rc != 0:
        raise ValueError(f"CM file {cm_path!r} failed cmstat's format check "
                         f"(rc={rc}): {stderr.strip()}")


def check_cm_source_formats(source):
    """check_cm_format for a --canonical-cm/--cyto-cm-db source: a single CM
    file, or a directory (checked recursively via find_cm_files)."""
    paths = find_cm_files(source) if os.path.isdir(source) else [source]
    for path in paths:
        check_cm_format(path)


def find_cm_files(cm_dir):
    """recursively list all .cm files under cm_dir."""
    return [
        os.path.join(root, f)
        for root, _, files in os.walk(cm_dir)
        for f in files if f.endswith(".cm")
    ]


def _scan_cm_files(cm_dir, pattern, key_fn, kind, exclude=None, warn_on_conflict=False):
    """shared walk-and-regex-match skeleton for CM index builders. key_fn(match)
    turns a regex match into the index key; files that don't match (or match
    `exclude`) are skipped with a debug log rather than guessed at, since
    mis-binning here would silently route to the wrong model. warn_on_conflict
    logs when a later file overwrites an earlier one under the same key, since
    that silently prefers one file over another rather than erroring."""
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


# --- cmalign: one call per (sequence, CM) ---

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
    because arm-loss diagnosis and CM-selection scoring need alignment-column
    coordinates."""
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
    (e.g. arm loss); nulling only that side and stripping would let naive
    re-matching silently re-pair the orphan with an unrelated stem (observed
    corrupting the acceptor stem). fix: read pairing from the full consensus db
    via RNA.ptable before stripping, and null BOTH sides of any pair where
    either column is gapped here, so stripping can't create an orphan at all."""
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


# --- stem/loop topology (forgi-based) ---

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

    Shared by get_stem_loop_elements (arm-loss diagnosis, mito-only) and
    parse_topology (Sprinzl labeling, shared), so both use the same physical
    stems."""
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


def find_anticodon_stem_index(aligned_seq, stem_loop_elements, anticodon, expected_index=None):
    """search for anticodon within each stem-loop's hairpin-loop columns only.
    both '-' and '.' must be stripped from loop sequences together; filtering
    only '-' can leave insert-column junk that produces a spurious extra match,
    breaking the "exactly one loop" assumption.

    expected_index breaks a tie between >=2 content matches, but only with
    >=3 stem-loops total: position alone identifies the arm there (D-arm
    first, anticodon-arm second), regardless of coincidental sequence
    matches elsewhere. with exactly 2 stem-loops it can't, since either one
    could be the arm that's missing.

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
    if (expected_index is not None and len(stem_loop_elements) >= 3
            and expected_index in candidates):
        return expected_index, f"positional_tiebreak_ambiguous_{len(candidates)}_loop_matches"
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


def _pick_by_anticodon_anchor(header, seq, anticodon, candidates):
    """given >=2 CM paths that could all serve the same (ambiguous) aa code,
    align to each and keep whichever anchors the header's anticodon in its
    own anticodon loop. the anticodon is the discriminating fact, not
    filename or dict-key suffix.

    - shared by mito's resolve_armless_cm (Leu1/Leu2, Ser1/Ser2 filenames) and
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
      case needs missing_arm, established via arm-loss diagnosis on the
      canonical alignment (mito.classify_arm_loss) when this sequence came
      from the mito path; always None on the cytosolic path, where 2-stem
      structures don't occur.
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
        loop_bases = "".join(seq[p] for p in direct_loop(idx, inner_stems[idx]) if seq[p] not in "-.")
        return ac and ac in loop_bases.upper().replace("T", "U")

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
        # min(after) breaks for class-ii tRNAs (ser, leu) and some tRNAs
        # with a variable arm stem: it picks the variable arm as t-arm instead.
        t_stem = max(after, key=lambda g: g["stem5_cols"][0]) if after else None
        # a real variable-ARM stem (class-ii: Leu, Ser) is whatever's left
        # in `after` besides t_stem; only trusted when exactly one such
        # candidate remains; with a single "after" candidate there's no way
        # to tell a bare variable arm from a missing T-arm from topology
        # alone, so that ambiguous case is left to the existing missing_arm
        # machinery (mito.classify_arm_loss) rather than guessed at here.
        v_candidates = [g for g in after if g is not t_stem]
        v_stem = v_candidates[0] if len(v_candidates) == 1 else None
        # a sequence that deleted every base of the arm has no arm to number
        if v_stem and not _occupied_count(seq, v_stem["stem5_cols"]
                                          + v_stem["loop_cols"] + v_stem["stem3_cols"]):
            v_stem = None

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
    - missing_arm (from mito.classify_arm_loss's diagnosis on the canonical
      alignment, if this sequence was rerouted) is passed through to
      locate_anticodon_stem to resolve its 2-stem case; see its docstring.
      always None on the cytosolic path."""
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
    such bulge in real tRNA data is a cmalign insert column (#=GC RF == '.'),
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


# --- alignment-column-aware Sprinzl assignment ---
#
# A CM's ss_cons marks every alignment column as match-state or insert-state,
# fixed for that CM regardless of which sequence is aligned to it: literal '.'
# is an insert column, anything else is a match column. sprinzl_map_from_alignment
# reads that per-column status straight off cmalign's raw output and assigns
# each block by walking its columns in 5'->3' order against two label pools:
# a core pool (the block's plain Sprinzl numbers) advanced by every
# match-state column, and an optional insertion pool (named codes like 17a)
# drawn on by occupied insert-state columns.
#
# Worked example, D-loop core positions 14-21 with insertion pools
# {"17": ["17a"], "20": ["20a", "20b"]}: four occupied match columns advance
# the core pool to 14, 15, 16, 17 in turn, each getting that label. An
# occupied insert column right after draws "17a" from 17's pool (the pool
# key is whichever core label the pointer last advanced to - not a running
# count of insertions seen so far). Two match-state deletions then advance
# the core pool through 18 and 19 silently: the pool still moves, but a
# deletion has no base to attach a label to. Two more occupied match columns
# take 20 and 21; an insert column right after 20 draws "20a" (
# the pointer has since moved past 20). Once both the core pool (8 slots)
# and the current label's insertion pool are exhausted, any further occupied
# column gets a generic letter suffix on the last label used (e.g. "21A").
#
# This all stays in raw (unstripped) column space, so parse_topology and
# locate_anticodon_stem take the CM's own ss_cons dot-bracket collapse
# directly: it is well-formed on its own, fixed per CM, with no per-sequence
# gap bookkeeping required to call them.


def _next_suffix(anchor, counts):
    """anchor + next unused letter (A, B, ..., Z, AA, ...); counts is shared
    across a block's calls so overflow letters for one anchor stay sequential."""
    n = counts.get(anchor, 0)
    counts[anchor] = n + 1
    letter = chr(ord("A") + n) if n < 26 else f"A{chr(ord('A') + n - 26)}"
    return f"{anchor}{letter}"


def _iter_block_labels(cols, core_slots, is_match, is_occupied, insertion_pools,
                        suffix_counts, anchor=None):
    """yield (col, label) for every occupied column in cols, running the core-pool
    / insertion-pool / overflow state machine described above.

    is_match(col) and is_occupied(col) classify each column; insertion_pools
    is {core_label: [reserved_codes]}; suffix_counts accumulates overflow
    letters per anchor label so repeated calls for one block stay sequential.

    anchor seeds the running label with the preceding block's last one, so an
    insert-state column at this block's 5' edge (reached before any core slot
    is consumed) still has something to suffix onto: a D-loop leading insert
    becomes 13A, off the last D-stem-5' position."""
    exhausted = object()
    core_iter = iter(core_slots)
    pool_used = defaultdict(int)
    # current stays the last CORE label, so consecutive overflows read 60A,
    # 60B, 60C rather than compounding into 60A, 60AA, 60AAA.
    current, label = anchor, None
    for col in cols:
        occ = is_occupied(col)
        if is_match(col):
            nxt = next(core_iter, exhausted)
            if nxt is not exhausted:
                current = nxt
                label = current if occ else None
            elif occ and current is not None:
                label = _next_suffix(current, suffix_counts)
            else:
                label = None
        elif occ and current is not None:
            pool = insertion_pools.get(current, ())
            used = pool_used[current]
            if used < len(pool):
                label, pool_used[current] = pool[used], used + 1
            else:
                label = _next_suffix(current, suffix_counts)
        else:
            label = None
        if label is not None:
            yield col, label


def _is_occupied(aligned_seq, col):
    return aligned_seq[col] not in "-."


def _occupied_count(aligned_seq, cols):
    return sum(1 for c in cols if _is_occupied(aligned_seq, c))


def _occupied_bases(aligned_seq, cols):
    return [aligned_seq[c].upper() for c in cols if _is_occupied(aligned_seq, c)]


def _raw_to_final_index(aligned_seq):
    """{raw_col: index_in_finalize_structure's_ungapped_seq} for every occupied
    column, in order - the coordinate translation assign_block needs to write
    labels keyed the same way sprinzl_map's caller-visible labels dict is."""
    raw_to_final, n = {}, 0
    for col, ch in enumerate(aligned_seq):
        if ch not in "-.":
            raw_to_final[col] = n
            n += 1
    return raw_to_final


def _split_occupied(aligned_seq, cols, n_head, n_tail):
    """split cols into (head, middle, tail): head/tail hold exactly n_head/
    n_tail occupied columns each; an interleaved gap column rides along with
    whichever side owns the adjacent occupied column."""
    occ_idx = [i for i, c in enumerate(cols) if _is_occupied(aligned_seq, c)]
    head_end = occ_idx[n_head - 1] + 1 if n_head else 0
    tail_start = occ_idx[-n_tail] if n_tail else len(cols)
    return cols[:head_end], cols[head_end:tail_start], cols[tail_start:]


def _assign_block(labels, cols, core_slots, ss_cons, aligned_seq, raw_to_final,
                   insertion_pools=None, suffix_counts=None, anchor=None):
    """write labels[raw_to_final[col]] for one Sprinzl block, per
    _iter_block_labels, reading match/insert state from ss_cons. returns the
    last label written, for the next block to anchor on (or anchor unchanged
    if this block wrote nothing)."""
    is_match = lambda col: ss_cons[col] != "."
    last = anchor
    for col, label in _iter_block_labels(
            cols, core_slots, is_match, lambda col: _is_occupied(aligned_seq, col),
            insertion_pools or {}, suffix_counts if suffix_counts is not None else {},
            anchor=anchor):
        labels[raw_to_final[col]] = label
        last = label
    return last


def _assign_plain_zip(labels, cols, core_slots, aligned_seq, raw_to_final,
                       suffix_counts=None, anchor=None):
    """assign core_slots, in order, to cols' occupied columns only, skipping
    gaps entirely rather than treating them as match-state deletions that
    advance the pool with no output. for blocks where match/insert-state
    carries no reliable Sprinzl meaning: the CCA trailer, which sits past the
    model's own consensus structure, and any block holding exactly as many
    bases as slots. returns the last label written, same contract as
    _assign_block."""
    suffix_counts = {} if suffix_counts is None else suffix_counts
    exhausted = object()
    core_iter = iter(core_slots)
    current = anchor
    last = anchor
    for col in cols:
        if aligned_seq[col] in "-.":
            continue
        nxt = next(core_iter, exhausted)
        if nxt is not exhausted:
            current = label = nxt
        elif current is not None:
            label = _next_suffix(current, suffix_counts)
        else:
            continue
        labels[raw_to_final[col]] = label
        last = label
    return last


def _assign_anticodon_loop_block(labels, cols, ss_cons, aligned_seq, raw_to_final,
                                  anticodon, suffix_counts, anchor=None):
    """locate the anticodon among the loop's occupied columns (centered-match,
    same anchoring as _assign_anticodon_loop), then run the flanking
    raw-column ranges through _assign_block so an indel there is handled by
    the same core/insertion-pool rule as every other block. returns the last
    label written, same contract as _assign_block."""
    ac = (anticodon or "").upper().replace("T", "U")
    occ_idx = [i for i, c in enumerate(cols) if _is_occupied(aligned_seq, c)]
    loop_seq = "".join(aligned_seq[cols[i]] for i in occ_idx).upper().replace("T", "U")
    matches = [m.start() for m in re.finditer(f"(?={re.escape(ac)})", loop_seq)] if ac else []
    if not matches:
        return _assign_block(labels, cols, [str(i) for i in range(32, 39)], ss_cons,
                             aligned_seq, raw_to_final, suffix_counts=suffix_counts,
                             anchor=anchor)

    center = (len(loop_seq) - 3) / 2
    ac_start = min(matches, key=lambda i: abs(i - center))
    ac_occ = occ_idx[ac_start:ac_start + 3]
    before, ac_cols, after = cols[:ac_occ[0]], cols[ac_occ[0]:ac_occ[-1] + 1], cols[ac_occ[-1] + 1:]

    n_before = _occupied_count(aligned_seq, before)
    before_slots = ([str(i) for i in range(34 - n_before, 34)] if n_before <= 2 else ["32", "33"])
    last = _assign_block(labels, before, before_slots, ss_cons, aligned_seq, raw_to_final,
                         suffix_counts=suffix_counts, anchor=anchor)
    for i, c in enumerate(c for c in ac_cols if _is_occupied(aligned_seq, c)):
        if i < 3:
            labels[raw_to_final[c]] = str(34 + i)
            last = str(34 + i)
    return _assign_block(labels, after, ["37", "38"], ss_cons, aligned_seq, raw_to_final,
                         suffix_counts=suffix_counts, anchor=last)


# a D-loop shorter than its eight slots loses length at the dihydrouridine
# positions, so 16, 17 and 20 give their bases up before 14, 15, 18, 19 and 21,
# which Biela et al. 2023 lists among the nucleotides conserved across tRNAs.
# Order taken from Suzuki et al. 2020's curated human mt-tRNAs.
D_LOOP_DROP_ORDER = ["17", "20", "16", "19", "18"]

# the D-loop takes its extra bases at 20a/20b far more often than at 17a
D_LOOP_INSERTION_ORDER = ["20a", "20b", "17a"]

# a short T-loop empties from its ends inward, alternating 3' then 5', and keeps
# the 56-58 core. Reproduces every T-loop in Suzuki et al. 2020's curated human
# mt-tRNAs.
T_LOOP_DROP_ORDER = ["60", "54", "59", "55"]

# 44, 45, 46 and 48 each hold a tertiary contact: G26-A44, G10-C25-G45,
# C13-G22-G46 and the Levitt pair G15-C48 (Biela et al. 2023). 47 holds none and
# gives its base up first. Nothing orders the other four against each other, so
# a block needing a second drop is left to the model.
V_LOOP_DROP_ORDER = ["47"]

SLOT_DROP_ORDER = {"d_loop": D_LOOP_DROP_ORDER, "t_loop": T_LOOP_DROP_ORDER,
                   "v_loop": V_LOOP_DROP_ORDER}

# slots either side of the conserved 18-19 pair, spent outward from it
D_LOOP_5P_SLOTS = ["14", "15", "16", "17", "17a"]
D_LOOP_3P_SLOTS = ["20", "20a", "20b"]


def _d_loop_slots_from_gg(bases):
    """slots for a D-loop, seating the GG that Biela et al. 2023 report at 18
    and 19 on those two positions and spending the rest outward from it. None
    when no GG sits close enough to 18 to be that pair.

    Counting slots from 14 instead gets the common eight-base loop wrong: it
    fills 14-21 solid, which lands the GG on 17-18, because 17 stays empty
    until a loop is long enough to need it."""
    n = len(bases)
    for k in range(3, min(len(D_LOOP_5P_SLOTS), n - 2) + 1):
        if bases[k] != "G" or bases[k + 1] != "G":
            continue
        after = n - k - 2
        if after > len(D_LOOP_3P_SLOTS) + 1:
            return None
        # 21 closes the loop against the D-stem, so it takes the last base and
        # 20/20a/20b fill the gap left in front of it
        tail = D_LOOP_3P_SLOTS[:after - 1] + ["21"] if after else []
        return D_LOOP_5P_SLOTS[:k] + ["18", "19"] + tail
    return None


def _shrink_slots(core_slots, n_bases, drop_order):
    """core_slots cut down to n_bases entries, dropping in drop_order and
    keeping the rest in Sprinzl order. Leaving the choice to the CM instead puts
    the gap wherever its deletions fell, which strands a conserved position."""
    # an order may name slots this block does not have, and those must not
    # count toward the number dropped
    candidates = [slot for slot in drop_order if slot in core_slots]
    dropped = set(candidates[:len(core_slots) - n_bases])
    return [slot for slot in core_slots if slot not in dropped]


def _expand_slots(core_slots, n_bases, insertion_pools, code_order=()):
    """core_slots grown toward n_bases by spending reserved insertion codes at
    their anchors. Bases beyond the core count give the number of insertions,
    so every core slot keeps a base and the spares take the reserved codes.
    code_order picks which codes go first; anything left over follows in
    anchor order."""
    extra = n_bases - len(core_slots)
    available = [(slot, code) for slot in core_slots for code in insertion_pools.get(slot, ())]
    ranked = sorted(available, key=lambda sc: (code_order.index(sc[1])
                                               if sc[1] in code_order else len(code_order)))
    spend = {code for _, code in ranked[:max(extra, 0)]}
    grown = []
    for slot in core_slots:
        grown.append(slot)
        grown.extend(code for code in insertion_pools.get(slot, ()) if code in spend)
    return grown


def _absorb_unclaimed_columns(specs):
    """extend each block to cover every column up to the next block's start,
    given specs already sorted by start column.

    forgi reports a stem's paired columns only, so a stem-internal bulge lands
    in no block's own column list and would otherwise go unlabeled. Handing it
    to the block it sits inside puts it through the same insertion rule as any
    other unpaired column, which suffixes it onto the label before it."""
    out = []
    for i, (cols, core_slots, pools, mode) in enumerate(specs):
        if i + 1 < len(specs):
            cols = list(range(cols[0], specs[i + 1][0][0]))
        out.append((cols, core_slots, pools, mode))
    return out


def sprinzl_map_from_alignment(alignment, anticodon, missing_arm=None, wc=False, header=""):
    """assign a Sprinzl label to every occupied column by reading match/
    insert/deletion status directly off cmalign's raw output (see module
    note above).

    - alignment: cmalign_one's return dict (raw aligned_seq/ss_cons, gapped).
    - anticodon, missing_arm: same meaning as sprinzl_map.
    - wc: how far a stem may be re-seated by base-pairing first (see
      slide_stems_to_improve_pairing); 0 skips it. sliding is by occupied
      columns, so a step moves the helix one base.
    - a sequence whose structure did not come from cmalign has no match/insert
      state to read and belongs on sprinzl_map instead; see mito's
      threading-failure branch.
    - returns {final_seq_index: label}, where final_seq_index matches
      finalize_structure's ungapped/uppercased seq (same base order) - pair
      with finalize_structure(alignment) for final_seq/final_ss."""
    aligned_seq, ss_cons = alignment["aligned_seq"], alignment["ss_cons"]
    raw_to_final = _raw_to_final_index(aligned_seq)
    ss_cons = slide_stems_in_alignment(aligned_seq, ss_cons, max_slide=wc, header=header)
    raw_db = drop_orphan_brackets(RNA.db_from_WUSS(ss_cons))
    topo = parse_topology(raw_db)
    arms = locate_anticodon_stem(topo, raw_db, aligned_seq, anticodon, missing_arm)

    specs = []

    def block(cols, core_slots, insertion_pools=None, mode=None):
        if cols:
            specs.append((list(cols), core_slots, insertion_pools, mode))

    def zip_block(cols, core_slots):
        if cols:
            specs.append((list(cols), core_slots, None, "zip"))

    block(topo["acceptor_5"], [str(i) for i in range(1, 8)])
    block(topo["acceptor_3"], [str(i) for i in range(66, 73)])
    zip_block(topo["trailer"], ["73", "74", "75", "76"])

    d_loop_pools = {"17": ["17a"], "20": ["20a", "20b"]}
    if arms["d_stem5"]:
        block(arms["linker_5"], ["8", "9"])
        block(arms["d_stem5"], [str(i) for i in range(10, 14)])
        block(arms["d_loop"], [str(i) for i in range(14, 22)], d_loop_pools, mode="d_loop")
        block(arms["d_stem3"], [str(i) for i in range(22, 26)])
        block(arms["linker_dc"], ["26"])
    else:
        # d-armless: the replacement loop occupies positions 8-26 in one run,
        # with the same reserved D-loop insertion codes at the same anchors.
        block(arms["linker_5"], [
            "8", "9", "10", "11", "12", "13",
            "14", "15", "16", "17", "18", "19", "20", "21",
            "22", "23", "24", "25", "26",
        ], d_loop_pools)

    block(arms["c_stem5"], [str(i) for i in range(27, 32)])
    if arms["c_loop"]:
        specs.append((list(arms["c_loop"]), None, None, "anticodon"))
    block(arms["c_stem3"], [str(i) for i in range(39, 44)])

    if arms["v_stem5"]:
        block(arms["ct_linker"], ["44", "45"], mode="v_loop")
        block(arms["v_stem5"], [f"e1{i}" for i in range(1, 8)])
        block(arms["v_loop"], [f"e{i}" for i in range(1, 6)])
        n3 = min(_occupied_count(aligned_seq, arms["v_stem3"]), 7)
        block(arms["v_stem3"], [f"e2{k}" for k in range(n3, 0, -1)])
        block(arms["vt_linker"], ["46", "47", "48"], mode="v_loop")
    elif _occupied_count(aligned_seq, arms["var_loop"]) > 5:
        before, middle, after = _split_occupied(aligned_seq, arms["var_loop"], 2, 3)
        block(before, ["44", "45"], mode="v_loop")
        block(middle, [f"e{i}" for i in range(1, 6)])
        block(after, ["46", "47", "48"], mode="v_loop")
    else:
        block(arms["var_loop"], ["44", "45", "46", "47", "48"], mode="v_loop")

    block(arms["t_stem5"], [str(i) for i in range(49, 54)])
    block(arms["t_loop"], [str(i) for i in range(54, 61)], mode="t_loop")
    block(arms["t_stem3"], [str(i) for i in range(61, 66)])

    specs.sort(key=lambda spec: spec[0][0])
    specs = _absorb_unclaimed_columns(specs)

    labels, suffix_counts, anchor = {}, {}, None
    for cols, core_slots, pools, mode in specs:
        # a block holding exactly as many bases as it has slots has only one
        # consistent labelling, so the CM's view of which columns are
        # matches carries no extra information there - and acting on it does
        # harm when the CM threaded the block poorly, which mt-tRNA loops
        # frequently do (bases parked in insert columns while the consensus
        # columns beside them are called deletions). read match/insert state
        # only where the counts disagree and the placement is in question.
        n_bases = _occupied_count(aligned_seq, cols)
        if mode == "d_loop":
            anchored = _d_loop_slots_from_gg(_occupied_bases(aligned_seq, cols))
            if anchored:
                core_slots = anchored
        drop_order = SLOT_DROP_ORDER.get(mode)
        if drop_order and core_slots:
            if n_bases < len(core_slots):
                core_slots = _shrink_slots(core_slots, n_bases, drop_order)
            elif n_bases > len(core_slots):
                core_slots = _expand_slots(core_slots, n_bases, pools or {},
                                           D_LOOP_INSERTION_ORDER)
        exact_fit = core_slots is not None and n_bases == len(core_slots)
        if mode == "anticodon":
            anchor = _assign_anticodon_loop_block(labels, cols, ss_cons, aligned_seq,
                                                   raw_to_final, anticodon, suffix_counts,
                                                   anchor=anchor)
        elif mode == "zip" or exact_fit:
            anchor = _assign_plain_zip(labels, cols, core_slots, aligned_seq, raw_to_final,
                                        suffix_counts, anchor=anchor)
        else:
            anchor = _assign_block(labels, cols, core_slots, ss_cons, aligned_seq,
                                    raw_to_final, pools, suffix_counts, anchor=anchor)
    return labels
# --- stem register correction (--wc) ---

# a helix seated further off than this is a threading failure, handled by
# mito's RNAfold patch instead
MAX_STEM_SLIDE = 2


def slide_offsets(max_slide):
    """offsets to try, nearest first, so the smallest move that gains wins."""
    return sorted([d for d in range(-max_slide, max_slide + 1) if d],
                  key=lambda d: (abs(d), d))


def _stem_pairs(ss, group):
    """(5' col, 3' col) for each paired column of one stem group, read from the
    pair table so a merged stem with a bulge keeps its true partners."""
    pt = RNA.ptable(ss)
    return [(i, pt[i + 1] - 1) for i in group["stem5_cols"] if pt[i + 1] > i + 1]


def _pair_cols(pairs):
    return {col for pair in pairs for col in pair}


def _count_wc(seq, pairs):
    return sum(1 for i, j in pairs if (seq[i], seq[j]) in WC_PAIRS)


def wuss_stems(ss_cons):
    """internal stems as [{'pairs': [(5' col, 3' col), ...]}, ...], 5'->3'.
    e.g. [{'pairs': [(10,25), (11,24), (12,23)]}, {'pairs': [(49,65), (50,64), (51,63)]}]

    WUSS marks the acceptor stem '(' ')' and every internal stem '<' '>', so
    the arms come from the consensus line. 
    Nested pairs belong to one stem; a disjoint span
    starts the next."""
    stack, pairs = [], []
    for col, sym in enumerate(ss_cons):
        if sym == "<":
            stack.append(col)
        elif sym == ">" and stack:
            pairs.append((stack.pop(), col))
    stems = []
    for i, j in sorted(pairs):
        if stems and i < stems[-1]["pairs"][-1][1]:
            stems[-1]["pairs"].append((i, j))
        else:
            stems.append({"pairs": [(i, j)]})
    return stems


def _occupied_cols(aligned_seq):
    return [c for c, ch in enumerate(aligned_seq) if ch not in "-."]


def _column_offset_for_bases(aligned_seq, edge, steps, direction, taken):
    """columns from `edge` out to the `steps`-th free base beyond it, in
    `direction`.

    Only columns holding a base count, so deletions are stepped over, and
    only columns no other helix owns, so the count never runs through a
    neighbouring stem. Returns None when another helix or the end of the
    sequence arrives first: with no free base to move onto, there is nothing
    to slide."""
    seen, col = 0, edge + direction
    while 0 <= col < len(aligned_seq):
        if col in taken:
            return None
        if aligned_seq[col] not in "-.":
            seen += 1
            if seen == steps:
                return abs(col - edge)
        col += direction
    return None


def slide_stems_in_alignment(aligned_seq, ss_cons, max_slide=1, header=""):
    """re-seat internal stems on the consensus line; returns a new ss_cons.

    The whole helix moves by one column offset, so it stays a helix. The
    offset is measured in free bases rather than columns: one step lands on
    the next base that no other helix owns, stepping over deletions.

    The anticodon stem stays put, since the numbering is anchored to it, and
    the acceptor stem is '(' ')' in WUSS so it is never a candidate. A move
    needs a strict gain in WC/wobble pairs. max_slide of 0 disables sliding;
    ss_cons also comes back unchanged when the anticodon stem cannot be
    identified by stem count."""
    if max_slide < 1:
        return ss_cons
    stems = wuss_stems(ss_cons)
    frozen_idx = EXPECTED_ANTICODON_ARM_INDEX.get(len(stems))
    if frozen_idx is None:
        logger.debug(f"{header}: {len(stems)} internal stems, cannot tell which is the "
                     f"anticodon arm; leaving the structure alone")
        return ss_cons

    paired = {c for c, sym in enumerate(ss_cons) if sym in "<>()"}
    out = list(ss_cons)
    for idx, stem in enumerate(stems):
        if idx == frozen_idx:
            continue
        pairs = stem["pairs"]
        cols = {c for p in pairs for c in p}
        taken = paired - cols
        base = _count_wc(aligned_seq, pairs)
        for steps in range(1, max_slide + 1):
            moved = _best_slide(aligned_seq, pairs, min(cols), max(cols), steps, taken, base)
            if moved is None:
                continue
            for i, j in pairs:
                out[i] = out[j] = ","
            for i, j in moved:
                out[i], out[j] = "<", ">"
            logger.info(f"{header}: re-seated a stem by {steps} base(s) "
                        f"({base} -> {_count_wc(aligned_seq, moved)} of {len(pairs)} "
                        f"pairs complementary)")
            break
    return "".join(out)


def _best_slide(aligned_seq, pairs, edge5, edge3, steps, taken, base):
    """moved pairs for a `steps`-base slide either way, or None if neither
    direction is reachable or improves on `base`."""
    for direction, edge in ((-1, edge5), (1, edge3)):
        delta = _column_offset_for_bases(aligned_seq, edge, steps, direction, taken)
        if delta is None:
            continue
        moved = [(i + direction * delta, j + direction * delta) for i, j in pairs]
        cols = {c for p in moved for c in p}
        if cols & taken or min(cols) < 0 or max(cols) >= len(aligned_seq):
            continue
        if _count_wc(aligned_seq, moved) > base:
            return moved
    return None


def slide_stems_to_improve_pairing(seq, ss, anticodon, missing_arm=None, header="",
                                    max_slide=1):
    """re-seat stems that pair better a position or two along; returns the
    corrected dot-bracket structure.

    For a gap-free structure, where one position is one base. The alignment
    equivalent is slide_stems_in_alignment; mito's RNAfold-patched arms come
    here instead, having a fold but no alignment behind them.

    Both strands move by one offset, so a stem keeps its length, bulges and
    loop and only changes position. A slide needs a strict gain in WC/wobble
    pairs, takes the smallest offset that gives one, and may not land on
    another stem's columns. max_slide bounds how far it may travel. The
    acceptor and anticodon stems stay put, since the rest of the numbering is
    anchored to them.

    Human MT-TS1 motivates this: every canonical CM in the library seats its
    D-stem with two mismatched pairs, one base off a fully paired register."""
    topo = parse_topology(ss)
    arms = locate_anticodon_stem(topo, ss, seq, anticodon, missing_arm)
    frozen = set(topo["acceptor_5"] + topo["acceptor_3"] + arms["c_stem5"] + arms["c_stem3"])

    groups = _forgi_stem_groups(ss)
    owned = {c for g in groups for c in g["stem5_cols"] + g["stem3_cols"]}

    moves = []
    for group in groups:
        pairs = _stem_pairs(ss, group)
        cols = _pair_cols(pairs)
        if not pairs or cols & frozen:
            continue
        blocked = owned - cols
        for offset in slide_offsets(max_slide):
            moved = [(i + offset, j + offset) for i, j in pairs]
            new_cols = _pair_cols(moved)
            if min(new_cols) < 0 or max(new_cols) >= len(seq) or new_cols & blocked:
                continue
            if _count_wc(seq, moved) > _count_wc(seq, pairs):
                moves.append((pairs, moved))
                break

    out = list(ss)
    for pairs, moved in moves:
        logger.info(f"{header}: re-seated a stem by {moved[0][0] - pairs[0][0]:+d} "
                    f"({_count_wc(seq, pairs)} -> {_count_wc(seq, moved)} of "
                    f"{len(pairs)} pairs complementary)")
        for i, j in pairs:
            out[i] = out[j] = "."
    for _, moved in moves:
        for i, j in moved:
            out[i], out[j] = "(", ")"
    return "".join(out)
