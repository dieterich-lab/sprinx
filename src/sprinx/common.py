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


# --- stem register correction (--wc) ---

# a helix seated more than 2 positions off is a threading failure, handled by
# mito's RNAfold patch instead
MAX_STEM_SLIDE = 2

SLIDE_OFFSETS = sorted([d for d in range(-MAX_STEM_SLIDE, MAX_STEM_SLIDE + 1) if d],
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


def slide_stems_to_improve_pairing(seq, ss, anticodon, missing_arm=None, header=""):
    """re-seat stems that pair better one or two positions along; returns the
    corrected dot-bracket structure.

    Both strands move by one offset, so a stem keeps its length, bulges and
    loop and only changes position. A slide needs a strict gain in WC/wobble
    pairs, takes the smallest offset that gives one, and may not land on
    another stem's columns. The acceptor and anticodon stems stay put, since
    the rest of the numbering is anchored to them.

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
        for offset in SLIDE_OFFSETS:
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
