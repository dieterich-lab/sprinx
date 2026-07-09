#!/usr/bin/env python3
"""
sprinx.py -- Sprinzl-coordinate annotation for mt-tRNAs.

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
   this script avoids it by never comparing scores across models at all.

3. pipeline
   a. align to a canonical CM with cmalign --notrunc --nonbanded -g. --canonical-cm
      accepts multiple sources tried in priority order (e.g. bacterial whole-family
      CM, then a metazoan per-AA directory); the first tier whose anticodon anchors
      unambiguously and threads its own stem to the full canonical 5bp wins,
      but a short thread never disqualifies a tier outright -- a real anticodon
      stem can genuinely be shorter than 5bp, so the best-threaded anchored tier
      is kept and used if no later tier does better. never by score/E-value,
      see (2) -- a CM built for the wrong clade can fail to thread a divergent
      sequence at all, or mis-thread the otherwise-invariant anticodon stem
      specifically even when the anticodon itself anchors. details in
      select_cm_and_align.
   b. anchor on the anticodon; a missing UPSTREAM arm (D-arm) shifts remaining
      structure into wrong model columns (register shift). missing DOWNSTREAM
      arm (T-arm) does not shift. measure offset = expected_anticodon_slot - observed.
   c. n_pairs==0 at a stem slot means zero alignment columns have BOTH pairing
      partners simultaneously non-gap. no column can form a pair, so no stem can
      exist there -- geometry forces the call, with no threshold to tune.
      n_pairs==0 has two distinct causes that require different responses:
        (i)  genuine arm loss: the sequence simply has no arm. the element span
             across the alignment is mostly or entirely gap characters.
        (ii) CM threading failure: the arm exists but cmalign placed its sequence
             into unmodeled insert columns because the arm is too divergent from
             the CM consensus. the stem model columns are all gaps, but the span
             DOES contain nucleotides as insert characters.
      distinguishing (i) from (ii): count non-gap nucleotides across the full
      element span (stem + loop model columns + intervening insert characters).
      if the count is < n_stem_cols + MIN_HAIRPIN_LOOP (=3, steric minimum for
      the RNA backbone to close a hairpin), no hairpin can form physically:
      genuine arm loss (i). otherwise: threading failure (ii).
      hybrid Infernal + RNAfold design: for threading failures, Infernal's
      canonical CM is correct for all arms it DID thread properly; only the
      mis-threaded arm needs structural recovery. RNAfold MFE on the short arm
      span (typically 13-20 nt) is reliable at this length because competing folds
      are energetically negligible. the hybrid avoids two failure modes: (a)
      relying on Infernal alone would call threading failures as arm loss and
      misroute to an armless CM; (b) relying on RNAfold alone for full-sequence
      mt-tRNA folding is unreliable due to tertiary interactions and base
      modifications not captured by 2D MFE.
   d. if genuinely absent, reroute to armless CM (Ozerova et al. 2024).
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
   is unaffected by which canonical CM tier won above -- it only triggers once
   a genuine arm-loss diagnosis is made from whichever tier's alignment was used.
   each --canonical-cm source is a directory of {label}_{AA}.cm files (e.g.
   Metazoan_P.cm; label/clade is ignored, selection is by AA only, per-sequence,
   same as armless CM selection) or a single CM file (applies to every aa,
   e.g. a whole-family CM like TRNAinf-bact.cm).

5. output
   sprinzl_mapping.tsv: seq_id, seq_index, nucleotide, sprinzl_position, region,
   cm_used, rerouted, arm_loss_call. optional --plot for an R2DT-rendered 2D
   diagram (one panel per sequence, stitched into a single SVG); see README
   for the R2DT Singularity image setup.

usage: see README.md, or `python sprinx.py --help`.
"""

import argparse
import glob
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import warnings
import xml.etree.ElementTree as ET
from collections import defaultdict

import cairosvg
import pandas as pd
import RNA
from forgi.graph.bulge_graph import BulgeGraph
from Bio.Data.IUPACData import protein_letters_3to1
from Bio import SeqIO
from loguru import logger
from scipy.stats import binomtest

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
# topological fact, not tunable -- changing it requires a different CM.
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
    None; returns None rather than guessing on format mismatch -- a wrong
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
    are skipped with a debug log rather than guessed at -- mis-binning here would
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
    label (e.g. clade) is ignored -- selection is by AA only. armless CM files
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
    (e.g. genuine arm loss); nulling only that side and stripping would let
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
    """all physical stems via forgi BulgeGraph, merging pieces a single interior-
    loop bulge splits apart (forgi has no notion of "one bulged helix" vs "two
    helices", so that merging has to happen here).
    an interior-loop (forgi 'i') edge connects exactly two stems, and its graph
    shape is the same whether it is a real bulge inside one helix (merge) or a
    junction between two distinct helices (don't merge). the junction case only
    arises when a multi-branch loop degenerates to two stems -- an armless CM
    leaving just the acceptor + C-stem (doubly-armless), or C-stem + T-stem; a
    full cloverleaf joins its arms through a multiloop ('m'), never 'i', so its
    arms never get merged. the tell is how many arm/acceptor "anchors" a merged
    group holds: each hairpin is one arm, and the outermost stem (lowest start
    column = acceptor) counts as one more anchor when the group also holds a
    hairpin. more than one anchor means the group spans a real junction, so keep
    the stems separate; one or none is a single-helix bulge (acceptor-internal
    or arm-internal), so merge. TestForgiStemGroups covers the four classes.
    each dict: {'stem5_cols', 'stem3_cols', 'stem_cols' (both sides), 'loop_cols'
    (hairpin loop, or [] for the acceptor / any hairpin-less stem), 'span'},
    0-indexed, sorted by span start. shared by get_stem_loop_elements (arm-loss
    diagnosis) and parse_topology (Sprinzl labeling) so both see the same
    physical stems instead of two stem parsers drifting apart."""
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
    G-U wobble pairs. p_value: binomial test vs null rate len(WC_PAIRS)/16, not
    thresholded (short D-arm stems often lack power even when real) -- callers
    read per_stem_complementarity directly instead of a binary verdict.
    raw WUSS in ss is handled transparently by db_from_WUSS."""
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
    p = binomtest(k, n, p=len(WC_PAIRS) / 16, alternative="greater").pvalue if n > 0 else None
    return {"n_pairs": n, "n_compatible": k, "p_value": p}


def classify_arm_loss(header, aligned_seq, ss_cons,
                      expected_anticodon_index=EXPECTED_ANTICODON_STEM_INDEX):
    """top-level structural diagnosis for one cmalign'd sequence: which arm is
    missing, measured via register shift (D-arm, when it occurs) or per-slot
    absent() (D-arm when no shift occurs, T-arm always, via MIN_STEM_PAIRS --
    see its docstring for why this is a soft signal). always returns full
    diagnostics for every stem, even on ambiguous input. see TestCanonical36,
    TestTArmless, TestDArmless, TestBothArmlessMature for end-to-end
    validation against real alignments."""
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
        return per_stem[i]["n_pairs"] < MIN_STEM_PAIRS

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
        # variable-arm stem some CMs model) are reported but never load-bearing
        # on their own -- no armless CM exists to reroute a variable-arm loss to.
        other_missing = [i for i in range(idx + 1, t_arm_idx) if absent(i)]

        if d_absent and t_absent:
            # doubly-armless tRNAs show offset==0 because with both arms absent
            # the anticodon arm still lands in the expected model columns (no
            # single-arm register shift occurs).
            result["call"] = f"BOTH_ARMS_MISSING_slots={[d_arm_idx, t_arm_idx]}"
            result["missing_arm"] = "d_and_t"
        elif d_absent:
            # D-arm absent but NO register shift -- cmalign left the D-arm's own
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
# this is a judgment call, not a certainty -- it only flags *candidates* for
# arm loss; every candidate still has to pass the hard arm_span_has_enough_sequence
# check before any reroute happens.
MIN_STEM_PAIRS = 3

# full canonical anticodon-stem length. NOT a hard minimum -- a real anticodon
# stem can genuinely thread as few as ~3 pairs, so a short count is never
# grounds to reject a tier (that would misfire on real biology). used only as
# a preference signal in select_cm_and_align: among tiers that anchor the
# anticodon cleanly, keep checking for one that reaches this full count rather
# than settling for the first anchor, since a poorly-fitting CM can mis-thread
# even an always-present stem while still anchoring the anticodon itself.
# confirmed on real data (mt-Cys): TRNAinf-bact.cm threads only 3 of 5 pairs;
# Metazoa_C.cm threads all 5 for the identical sequence.
ANTICODON_STEM_PAIRS = 5


def arm_span_has_enough_sequence(aligned_seq, elem):
    """first-stage (fast, hard) filter after a stem slot is flagged absent: does
    the span contain enough nucleotides to physically form a hairpin
    (n_stem_cols + MIN_HAIRPIN_LOOP, the steric minimum)? False here means
    definite genuine loss. True is not proof of a real arm, just not ruled out by
    volume alone -- see arm_is_threading_failure for the required 2nd check."""
    start, end = elem["span"]
    n_nts = sum(1 for c in aligned_seq[start:end] if c not in "-.")
    return n_nts >= len(elem["stem_cols"]) + MIN_HAIRPIN_LOOP


def _arm_full_span_subseq_and_fold(aligned_seq, final_seq, elem):
    """extract the FULL non-gap span (matched + insert columns together) and
    fold it with RNAfold MFE. full span, not insert-only: a real arm's
    sequence can land in the slot's own matched columns too (e.g. human
    mt-Val's T-arm under TRNAinf-bact.cm), which are safe to fold over since
    this slot was already flagged absent (n_pairs below MIN_STEM_PAIRS) --
    no real base pairs there to protect. shared by arm_is_threading_failure
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
    sequence alone. True: real, recoverable arm. False: genuine loss despite
    passing the count check."""
    _, arm_ss = _arm_full_span_subseq_and_fold(aligned_seq, final_seq, elem)
    return arm_ss is not None and "(" in arm_ss


def patch_threading_failure_arm(header, aligned_seq, final_seq, final_ss, elem):
    """recover arm structure for a confirmed CM threading failure (called only
    after arm_is_threading_failure): splice elem's own RNAfold fold into
    final_ss rather than refolding the whole molecule -- full-sequence RNAfold
    on a mt-tRNA is unreliable (tertiary contacts, modified bases), but a short
    isolated span (13-20nt) is fine.
    the span is cleared before the fold is written in: cmalign's own structure
    inside a flagged-as-failed span is exactly what's untrustworthy here (that's
    why arm_is_threading_failure fired), so a stray weak bracket that survived
    it (below MIN_STEM_PAIRS but still non-'.') must not be able to block the
    very patch meant to replace it -- observed on real data (mt-Cys under
    TRNAinf-bact.cm), where a single leftover pair blocked a 3bp D-stem
    RNAfold found that agreed with it and simply extended it. a bracket
    OUTSIDE the span whose partner falls INSIDE it would dangle once the span
    is cleared, so that partner is cleared too. both overrides are logged as
    ONE consolidated warning (not one line per position, and not debug-only):
    this is RNAfold overruling cmalign's own structural call, not merely
    filling in a gap, and it can in principle reach outside the span this
    function is actually meant to patch -- header identifies which sequence,
    since worker processes interleave in a multi-sequence run.
    returns patched final_ss, or the original on no fold / unbalanced result
    (safety net)."""
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
            f"{header}: RNAfold patch overrode cmalign's own structure -- "
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
    align to each and keep whichever anchors the header's anticodon in its own
    anticodon loop -- the anticodon is the discriminating fact, not filename or
    dict-key suffix. shared by resolve_armless_cm (Leu1/Leu2, Ser1/Ser2
    filenames) and _resolve_canonical_for_tier (same ambiguity, but from a
    bare GtRNAdb-style aa field with no isoacceptor digit at all). falls back
    to the first candidate, logged, if none anchor or no anticodon is known."""
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
    """top-level CM selection for one sequence.
    1. try each canonical CM tier in order (e.g. bacterial whole-family CM
       first -- mitochondria's endosymbiotic origin -- then a metazoan per-AA
       directory); never by score/E-value across tiers (module docstring
       section 2). a clean anticodon anchor is necessary but not sufficient:
       the anticodon stem's own length (ANTICODON_STEM_PAIRS) is a structural
       invariant a poorly-fitting CM can still mis-thread even when the
       anchor itself is clean (see mt-Cys: TRNAinf-bact.cm threads 3/5,
       Metazoa_C.cm threads 5/5, identical sequence) -- but a real anticodon
       stem can also genuinely be shorter than 5bp, so a short thread is never
       grounds to reject a tier outright. instead: keep the anchored tier
       with the fullest thread seen among tiers tried so far, stop as soon as
       one reaches the full canonical count, and accept the best available if
       none do -- preferring better evidence when the tiers on offer actually
       provide it, never disqualifying a real anchor. a tier that doesn't
       apply to this aa, or whose cmalign fails, is skipped; if none anchor
       at all, fall back to the first tier that aligned.
    2. no arm missing (or only variable arm): return the canonical alignment.
    3. D-arm (no register shift) or T-arm flagged absent: cross-check with
       arm_span_has_enough_sequence then arm_is_threading_failure before
       trusting it; passing both means patch via RNAfold instead of
       rerouting. D-arm via register shift (offset>0) skips this and is
       trusted directly -- its span holds sequence displaced by the shift
       itself, not the D-arm's own (or absent) content, so both checks would
       false-positive there.
    4. genuine loss: reroute via resolve_armless_cm (anticodon-disambiguated);
       no matching armless CM: warn and keep the canonical alignment.
    canonical_cm_tiers: list of tiers (path, or {aa_code: path} dict), or a
    bare path/dict wrapped as a single tier. returns dict: final_alignment,
    cm_used, diagnosis, rerouted, threading_failure_elem."""
    if isinstance(canonical_cm_tiers, (str, dict)):
        canonical_cm_tiers = [canonical_cm_tiers]

    canonical_alignment = canonical_cm = diagnosis = best_c_pairs = None
    first_alignment = first_cm = first_diag = None  # ultimate fallback: no tier ever anchors
    for tier in canonical_cm_tiers:
        path = _resolve_canonical_for_tier(header, seq, tier)
        if path is None:
            logger.info(f"{header}: skipping a canonical CM tier -- no CM for this amino acid there")
            continue
        aln = cmalign_one(header, seq, path)
        if aln is None:
            logger.info(f"{header}: moving to next canonical CM tier -- alignment failed")
            continue
        diag = classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
        if first_alignment is None:                    # ultimate fallback if nothing ever anchors
            first_alignment, first_cm, first_diag = aln, path, diag
        idx = diag["anticodon_stem_index"]
        if idx is None:
            logger.warning(
                f"{header}: anticodon did not anchor cleanly against {path}, "
                f"moving to next canonical CM tier\n"
                f"  aligned_seq={aln['aligned_seq']}\n"
                f"  ss_cons={aln['ss_cons']}"
            )
            continue

        # clean anchor is necessary but not sufficient: a CM that threads the
        # anticodon stem short of its full canonical length (ANTICODON_STEM_PAIRS)
        # may just be a poor fit for this sequence at that stem specifically (see
        # mt-Cys: TRNAinf-bact.cm threads 3/5, Metazoa_C.cm threads 5/5 for the
        # identical sequence) -- but a real, correctly-threaded anticodon stem can
        # also genuinely be shorter than 5bp, so a short count is NOT grounds to
        # reject a tier outright (that would misfire on real biology). instead:
        # keep the anchored tier with the fullest thread seen so far, keep
        # checking remaining tiers only while short of the canonical count, and
        # stop as soon as one reaches it -- never disqualifying an anchor, only
        # preferring a fuller one when the tiers on offer actually provide one.
        n_pairs = diag["per_stem_complementarity"][idx]["n_pairs"]
        if canonical_alignment is None or n_pairs > best_c_pairs:
            canonical_alignment, canonical_cm, diagnosis, best_c_pairs = aln, path, diag, n_pairs
        if n_pairs >= ANTICODON_STEM_PAIRS:
            break
        logger.warning(
            f"{header}: anticodon anchored against {path} but its own stem threaded only "
            f"{n_pairs}/{ANTICODON_STEM_PAIRS} pairs, checking remaining tiers for a fuller thread\n"
            f"  aligned_seq={aln['aligned_seq']}\n"
            f"  ss_cons={aln['ss_cons']}"
        )

    if canonical_alignment is None:
        canonical_alignment, canonical_cm, diagnosis = first_alignment, first_cm, first_diag
    if canonical_alignment is None:
        return _routing_result(None, None, None)

    missing_arm = diagnosis["missing_arm"]

    if missing_arm not in ("d", "t", "d_and_t"):
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    # step 3 (see docstring). D-arm via register shift is trusted directly: its
    # span contains non-D-arm sequence placed there by the shift itself, so
    # both checks would false-positive on genuinely D-armless sequences.
    # D-arm via no-shift (offset==0) doesn't have that problem and gets the
    # same cross-check as T-arm.
    if missing_arm in ("t", "d"):
        if missing_arm == "d" and diagnosis["register_offset"] != 0:
            elem = None  # register-shift D-arm: trusted directly, no cross-check
        else:
            elements = get_stem_loop_elements(canonical_alignment["ss_cons"])
            if missing_arm == "t":
                elem = elements[-1] if elements else None
            else:
                d_idx = diagnosis["anticodon_stem_index"] - 1
                elem = elements[d_idx] if 0 <= d_idx < len(elements) else None

        if elem and arm_span_has_enough_sequence(canonical_alignment["aligned_seq"], elem):
            final_seq, _ = finalize_structure(canonical_alignment)
            if arm_is_threading_failure(canonical_alignment["aligned_seq"], final_seq, elem):
                logger.info(
                    f"{header}: CM diagnosed {missing_arm}-arm missing against {canonical_cm} "
                    f"({diagnosis['call']}) but the span folds as a real hairpin "
                    f"(CM threading failure, not genuine arm loss) -- patching via RNAfold\n"
                    f"  aligned_seq={canonical_alignment['aligned_seq']}\n"
                    f"  ss_cons={canonical_alignment['ss_cons']}"
                )
                return _routing_result(canonical_alignment, canonical_cm, diagnosis,
                                        threading_failure_elem=elem)

    # genuine arm loss: reroute
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
    """acceptor stem (the one stem with no hairpin loop of its own) + inner
    stems, via _forgi_stem_groups -- not a hand-rolled contiguous-bracket scan,
    so a bulge or a finalize_structure-nulled pair inside a stem (e.g. the
    acceptor's own 5' or 3' half) doesn't get split off as a phantom extra
    stem and silently dropped from Sprinzl labeling. does not label D/C/T
    yet -- that needs the anticodon (see locate_anticodon_stem)."""
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


def locate_anticodon_stem(topo, ss, seq, anticodon, missing_arm=None):
    """identify C-stem (anticodon arm) by anticodon content, D-stem = sibling
    before C that does not enclose it, T-stem = sibling after it. inner_stems
    is a list of _forgi_stem_groups dicts (stem5_cols/stem3_cols/loop_cols
    already known from forgi -- no re-deriving them via pair_table here).
    innermost-stem-first search order (smallest span first) prevents a
    D-armless pseudostem (which structurally encloses C) from being matched
    before the real C-stem. 'does not enclose' check: a pseudostem that opens
    before C and closes after C must be excluded as D-arm candidate -- it IS
    the enclosing pseudostem.
    see TestSprinzlAssignment::test_d_armless_replacement_loop_gets_d_arm_labels.
    when the anticodon substring happens to occur in more than one stem's loop
    (observed on a real T-armless armless-CM alignment, two short remaining
    loops, coincidental match in both), don't guess by span size -- use
    missing_arm, already established via the far more reliable register-offset
    diagnosis on the canonical alignment (classify_arm_loss), to break the tie:
    D always precedes C, so if the T-arm is the one missing here (only D and C
    remain) the later candidate is C; if the D-arm is missing (only C and T
    remain) the earlier one is C. this positional rule is only valid when
    exactly TWO stems remain -- so it is gated on len(inner_stems)==2, NOT on
    missing_arm alone: a canonical 3-stem alignment kept as an RNAfold-patch
    fallback still carries missing_arm=d/t from the canonical diagnosis (e.g.
    human mt-Val), and applying a 2-stem positional shortcut to a 3-stem
    structure would mislabel it. falls back to the smallest-span candidate
    (previous behaviour) whenever the tie can't be safely broken this way."""
    inner_stems = topo["inner_stems"]

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

    ac = (anticodon or "").upper().replace("T", "U")
    search_order = sorted(range(len(inner_stems)),
                          key=lambda i: inner_stems[i]["span"][1] - inner_stems[i]["span"][0])
    matches = [idx for idx in search_order
               if ac and ac in "".join(seq[p] for p in direct_loop(idx, inner_stems[idx]))]

    c_idx = None
    if len(matches) == 2 and len(inner_stems) == 2 and missing_arm in ("d", "t"):
        by_position = sorted(matches, key=lambda i: inner_stems[i]["stem5_cols"][0])
        c_idx = by_position[-1] if missing_arm == "t" else by_position[0]
        logger.debug(f"anticodon {ac!r} matched both remaining stems; "
                     f"disambiguated via missing_arm={missing_arm!r} -> stem {c_idx}")
    elif len(matches) > 1 and len(inner_stems) == 3:
        # full D-C-T cloverleaf with the anticodon substring coincidentally in
        # more than one loop: the C-stem is the middle one by position. trust
        # that only when the middle stem is itself a match (its loop holds the
        # anticodon, as the real C-loop must); otherwise fall through.
        middle = sorted(range(3), key=lambda i: inner_stems[i]["stem5_cols"][0])[1]
        c_idx = middle if middle in matches else matches[0]
        logger.debug(f"anticodon {ac!r} matched {len(matches)} loops in a 3-stem "
                     f"structure; taking the middle stem -> {c_idx}")
    elif matches:
        c_idx = matches[0]
    elif inner_stems:
        c_idx = search_order[0]
    c_stem = inner_stems[c_idx] if c_idx is not None else None

    d_stem, t_stem, v_stem = None, None, None
    # outermost (last) column of the c-stem 3' strand -- stem3_cols[0] is the
    # INNERMOST column (adjacent to the loop); using it here was a real bug:
    # var_loop's boundary must start after the whole c-stem ends, not after
    # its first (innermost) column, or var_loop's own assign_slots call
    # silently overwrites the c-stem-3 columns between [0] and [-1] with wrong
    # (var-loop) labels -- confirmed on real data (e.g. mt-Glu): columns that
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
        # a genuine variable-ARM stem (class-ii: Leu, Ser) is whatever's left
        # in `after` besides t_stem -- only trusted when exactly one such
        # candidate remains; with a single "after" candidate there's no way
        # to tell a bare variable arm from a missing T-arm from topology
        # alone, so that ambiguous case is left to the existing missing_arm
        # machinery (classify_arm_loss) rather than guessed at here.
        v_candidates = [g for g in after if g is not t_stem]
        v_stem = v_candidates[0] if len(v_candidates) == 1 else None

    # outermost D-stem 3' column (strand edge), so the connector (pos 26) starts
    # only after the whole D-stem 3' strand -- a D-stem-internal 3' bulge sits
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
    # the plain 44-48 sequential run -- see sprinzl_map. ct_linker/vt_linker
    # are the (up to 2 / up to 3) genuinely unpaired nt flanking the v-stem
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
    shift the anticodon off 34-35-36 -- the tool's core deliverable -- even
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
    """assign a Sprinzl label to every nucleotide index; returns {seq_index: label}.
    D-armless tRNAs: replacement loop (all of linker_5) is mapped onto D-arm Sprinzl
    positions 8-26 by structural analogy, following Ozerova et al. 2024.
    missing T-arm produces no labels for its region. missing_arm (from
    classify_arm_loss's diagnosis on the canonical alignment, if this sequence was
    rerouted) is passed through to locate_anticodon_stem to break a same-content
    anticodon-match tie -- see its docstring."""
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
        # e-series -- see locate_anticodon_stem's v_stem docstring. reserved
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
    else:
        assign_slots(labels, arms["var_loop"], ["44", "45", "46", "47", "48"])
    assign_slots(labels, arms["t_stem5"],    [str(i) for i in range(49, 54)])
    assign_slots(labels, arms["t_loop"],     [str(i) for i in range(54, 61)])
    assign_slots(labels, arms["t_stem3"],    [str(i) for i in range(61, 66)])

    # strand ranges (first..last paired column of each stem strand) -- a
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
    construction) is a real nucleotide with no named Sprinzl slot -- letter-
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
    """worker for one (header, seq) FASTA record. takes a single tuple for
    Pool.map compatibility; canonical_cm_tiers and armless_cm_index are inside
    the tuple because each worker is a fresh process and module-level globals
    aren't reliably shared across fork vs spawn. per-tier canonical CM
    resolution (which .cm path applies to this aa, if any) happens inside
    select_cm_and_align, not here -- see _resolve_canonical_for_tier."""
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
    cm_only_ss = None
    if routing.get("threading_failure_elem"):
        cm_only_ss = final_ss
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
                     f"n_compatible={stem['n_compatible']} p={stem['p_value']}")
    if routing.get("threading_failure_elem"):
        logger.debug(f"  threading failure: patched via RNAfold; ss: {final_ss}")
    if routing["rerouted"]:
        logger.debug(f"  rerouted to: {routing['cm_used']}")
    logger.debug(f"  raw stockholm:\n{alignment['raw_sto']}")

    rows = []
    for i, base in enumerate(final_seq):
        label = sprinzl.get(i, "")
        # 'e'-prefixed variable-arm labels (e11, e1, e23, ...) aren't purely
        # numeric like every other slot -- match the optional 'e' along with
        # the digits so an overflow-suffixed one (e17A) still resolves to its
        # base code (e17) for the region lookup, same as '60A' -> '60' does.
        region_key = re.match(r"e?\d+", label).group() if label else ""
        rows.append({
            "seq_id": header, "seq_index": i, "nucleotide": base,
            "sprinzl_position": label, "region": SPRINZL_REGION.get(region_key, ""),
            "cm_used": cm_name, "rerouted": routing["rerouted"],
            "arm_loss_call": diagnosis.get("call"),
        })

    logger.info(f"{header}: {summary}  [{diagnosis.get('call')}]")
    return {"header": header, "rows": rows, "summary": summary,
            "seq": final_seq, "ss": final_ss, "sprinzl": sprinzl,
            "cm_only_ss": cm_only_ss}


# --- plotting (optional, --plot) ---
#
# renders each record's OWN final structure (final_ss -- sprinx's arm-loss/
# threading-failure diagnosis already baked in) via R2DT's template-free
# "stockholm" mode, not a structure R2DT would re-derive itself from its own
# template library (which could silently disagree with sprinx's own call).
# R2DT only accepts a real multi-sequence alignment as input, and these
# records aren't aligned to each other at all, so build_r2dt_stockholm fakes
# one: every record's sequence is concatenated end-to-end into a single row,
# with one #=GC structureID region marking each record's own column span --
# see https://docs.r2dt.bio for the annotation format.

R2DT_DEFAULT_IMAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib", "r2dt")


def _sanitize_r2dt_name(text, maxlen=40):
    """structureID/regionID names may not contain '|' or '.' (the annotation
    line's own delimiters); collapse anything else to '_' and truncate."""
    safe = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return safe[:maxlen] or "seq"


def _r2dt_id_line(segments):
    """build one #=GC structureID/regionID line: '|' marks each segment's first
    column, the segment's name fills the following columns, '.' pads the rest --
    the format R2DT's stockholm parser expects (see module comment above)."""
    total = sum(length for _, length in segments)
    line = ["."] * total
    pos = 0
    for name, length in segments:
        line[pos] = "|"
        for i, ch in enumerate(name[:length - 1]):
            line[pos + 1 + i] = ch
        pos += length
    line[-1] = "|"
    return "".join(line)


def build_r2dt_stockholm(plotted):
    """records (each with 'header', 'seq', 'ss') -> (Stockholm text, region
    names) with one structureID region per record, in the same order as
    `plotted`. regionID is the aa field, so R2DT's --color-by region groups
    isoacceptors under one colour. names are de-duplicated with a numeric
    suffix since sanitizing distinct headers can collide."""
    names, seen = [], {}
    for r in plotted:
        base = _sanitize_r2dt_name(r["header"])
        n = seen.get(base, 0)
        seen[base] = n + 1
        names.append(f"{base}_{n}" if n else base)

    concat_seq = "".join(r["seq"] for r in plotted)
    concat_ss = "".join(r["ss"] for r in plotted)
    structure_id = _r2dt_id_line([(name, len(r["seq"])) for name, r in zip(names, plotted)])
    region_id = _r2dt_id_line([(header_to_aa(r["header"]) or "na", len(r["seq"])) for r in plotted])

    sto_text = "\n".join([
        "# STOCKHOLM 1.0", "",
        f"seqs {concat_seq}", "",
        f"#=GC SS_cons     {concat_ss}",
        f"#=GC structureID {structure_id}",
        f"#=GC regionID    {region_id}",
        "//",
    ]) + "\n"
    return sto_text, names


def _cm_only_plot_path(path):
    """cloverleaves.svg -> cloverleaves_CMonly.svg, for the pre-RNAfold-patch
    comparison plot alongside the regular one."""
    root, ext = os.path.splitext(path)
    return f"{root}_CMonly{ext}"


SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def _flip_panel_north(panel, width, height):
    """mirror a panel vertically in place (flip y only, x untouched): R2R (the
    template-free layout engine R2DT uses here) always draws the acceptor
    stem at the bottom, with no orientation flag exposed to change that --
    confirmed consistent across every sequence/shape checked, so a single
    unconditional flip suffices. mirror y only -- a full 180-degree rotation
    would negate x too and swing every side arm from east to west; mirroring y
    alone moves the acceptor stem to the top while leaving east/west exactly
    where R2R put them.
    every <text> glyph gets its own counter-mirror about its own y: two
    y-mirrors about different lines compose into a pure translation, so the
    glyph itself stays upright while still landing at its correctly-mirrored
    position -- only the backbone/pairing geometry actually flips."""
    wrapper = ET.Element(f"{{{SVG_NS}}}g", {"transform": f"matrix(1 0 0 -1 0 {height})"})
    for child in list(panel):
        panel.remove(child)
        wrapper.append(child)
    panel.append(wrapper)
    for text in wrapper.iter(f"{{{SVG_NS}}}text"):
        ty = text.get("y", "0")
        text.set("transform", f"matrix(1 0 0 -1 0 {2 * float(ty)})")


NUCLEOTIDE_LETTERS = set("ACGUN")
SPRINZL_LABEL_STEP = 5


def _inject_sprinzl_labels(panel, sprinzl, label_step=SPRINZL_LABEL_STEP):
    """replace R2DT's own plain sequence-position numbering (1, 2, 3, ...,
    shown every 10th residue by default -- unrelated to Sprinzl coordinates)
    with sprinx's own Sprinzl labels. shown every label_step-th INTEGER
    position, but always for lettered insertions (17a, 20a, ...) since those
    don't follow a regular numeric cadence and would otherwise never appear.
    each nucleotide is its own top-level <g><title>i (...)</title><text>BASE
    </text></g> emitted by R2DT in strict 5'->3' order, so a running count of
    real base letters (skipping the synthetic 5'/3' end markers) lines up
    exactly with sprinzl's own 0-indexed final_seq positions."""
    for g in list(panel):
        text = g.find(f"{{{SVG_NS}}}text")
        line = g.find(f"{{{SVG_NS}}}line")
        cls = (text.get("class") if text is not None else None) or \
              (line.get("class") if line is not None else None) or ""
        if "numbering-label" in cls or "numbering-line" in cls:
            panel.remove(g)

    seq_idx = 0
    for g in list(panel):
        text = g.find(f"{{{SVG_NS}}}text")
        if text is None or (text.text or "").strip() in ("5'", "3'", ""):
            continue
        if (text.text or "").strip().upper() not in NUCLEOTIDE_LETTERS:
            continue
        label = sprinzl.get(seq_idx, "")
        seq_idx += 1
        show = label and (
            not label[:-1].isdigit()
            or int(re.match(r"\d+", label).group()) % label_step == 0
            or label == "1"
        )
        if not show:
            continue
        x, y = float(text.get("x")), float(text.get("y"))
        label_el = ET.SubElement(panel, f"{{{SVG_NS}}}text", {
            "x": str(x - 10), "y": str(y + 13),
            "class": "numbering-label",
        })
        label_el.text = label


CAPTION_FONT_SIZE = 11
CAPTION_LINE_HEIGHT = 14


def _wrap_caption(text, cell_w, max_lines=4):
    """text (header and summary, '\\n'-separated) -> wrapped lines that fit
    cell_w, each original line wrapped independently so the header and
    summary never run together into one blob. width is estimated from
    monospace glyph width since this is drawn as SVG <text>, not measured
    by a real layout engine."""
    chars_per_line = max(int(cell_w / (CAPTION_FONT_SIZE * 0.62)), 10)
    lines = [line for para in text.split("\n")
             for line in (textwrap.wrap(para, width=chars_per_line) or [""])]
    return lines[:max_lines]


def _grid_svg(panels, ncols, gap=20):
    """panels: list of (svg_root_element, width, height, caption, sprinzl).
    arranges them into a real grid (ncols per row), unlike R2DT's own
    --stitch (which only lays panels left-to-right in a single row --
    unusable past a handful of sequences, since the file just gets wider
    without bound). cell size is uniform (max panel width/height across all
    panels) so rows/columns stay aligned; each panel keeps its own native
    size, centered in its cell. captions are wrapped to the cell width (see
    _wrap_caption) and centered above the panel, not left-aligned to the
    cell -- a narrow panel in a wide cell would otherwise read as belonging
    to its neighbour. nested <svg> elements are SVG's own mechanism for
    embedding one diagram inside another at a given position/size -- no
    rasterization needed to compose them."""
    cell_w = max(w for _, w, _, _, _ in panels) + gap
    caption_lines = [_wrap_caption(caption, cell_w) for _, _, _, caption, _ in panels]
    caption_height = max(len(lines) for lines in caption_lines) * CAPTION_LINE_HEIGHT + 6
    cell_h = max(h for _, _, h, _, _ in panels) + gap + caption_height
    nrows = -(-len(panels) // ncols)

    root = ET.Element(f"{{{SVG_NS}}}svg", {
        "width": str(ncols * cell_w), "height": str(nrows * cell_h),
    })
    for i, ((panel, w, h, _, sprinzl), lines) in enumerate(zip(panels, caption_lines)):
        row, col = divmod(i, ncols)
        x = col * cell_w + (cell_w - w) / 2
        y = row * cell_h + caption_height
        cx = x + w / 2
        for li, line in enumerate(lines):
            text = ET.SubElement(root, f"{{{SVG_NS}}}text", {
                "x": str(cx), "y": str(row * cell_h + 12 + li * CAPTION_LINE_HEIGHT),
                "font-family": "monospace", "font-size": str(CAPTION_FONT_SIZE),
                "text-anchor": "middle",
            })
            text.text = line
        _inject_sprinzl_labels(panel, sprinzl)
        _flip_panel_north(panel, w, h)
        panel.set("x", str(x))
        panel.set("y", str(y))
        root.append(panel)
    return ET.ElementTree(root)


def make_plot(results, out_path, r2dt_image=R2DT_DEFAULT_IMAGE, ncols=6):
    """R2DT-rendered 2D diagram, one panel per successfully-processed record,
    arranged into our own grid (see _grid_svg -- R2DT's own stitching doesn't
    wrap into rows). plotted in header order: species (taxon field) first,
    then tRNA (aa field), so isoacceptors of the same species group together
    and species cluster in the figure. runs R2DT via its Singularity image
    (see README for setup); failures are logged and skipped, not raised,
    since --plot is a sanity-check convenience, not the tool's actual output."""
    plotted = [r for r in results if r["rows"]]
    if not plotted:
        logger.warning("nothing to plot -- every record failed upstream")
        return
    plotted.sort(key=lambda r: (header_to_taxon(r["header"]) or "",
                                header_to_aa(r["header"]) or "", r["header"]))

    sto_text, names = build_r2dt_stockholm(plotted)
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "sprinx_plot.sto"), "w") as f:
            f.write(sto_text)
        cmd = ["singularity", "exec", "-B", f"{tmpdir}:/rna/r2dt/temp", r2dt_image,
               "r2dt.py", "stockholm", "/rna/r2dt/temp/sprinx_plot.sto",
               "/rna/r2dt/temp/out", "--no-stitch"]
        _, stderr, rc = run(cmd)
        if rc != 0:
            logger.warning(f"R2DT plotting failed: {stderr.strip()[:500]}")
            return

        svg_dir = os.path.join(tmpdir, "out", "results", "svg")
        panels = []
        for name, r in zip(names, plotted):
            candidates = glob.glob(os.path.join(svg_dir, f"{name}_*.svg"))
            if not candidates:
                logger.warning(f"{r['header']}: R2DT produced no diagram for this sequence, skipping")
                continue
            root = ET.parse(candidates[0]).getroot()
            width, height = float(root.get("width")), float(root.get("height"))
            panels.append((root, width, height, f"{r['header']}\n{r['summary']}", r["sprinzl"]))

        if not panels:
            logger.warning("R2DT plotting produced no SVG output")
            return

        grid_path = os.path.join(tmpdir, "grid.svg")
        _grid_svg(panels, ncols).write(grid_path)

        ext = os.path.splitext(out_path)[1].lower()
        if ext == ".svg":
            shutil.copy(grid_path, out_path)
        else:
            _convert_svg(grid_path, out_path, ext)


_SVG_CONVERTERS = {".png": cairosvg.svg2png, ".pdf": cairosvg.svg2pdf}


# cairo's hard surface-size limit is ~32767px/side; a wide stitched plot (many
# sequences) can exceed that at scale=2.0, so the PNG scale is capped to keep
# the longer side under this, well clear of the real limit.
MAX_PNG_DIM = 16000


def _svg_intrinsic_size(svg_path):
    """(width, height) in px from an SVG's own width/height attributes."""
    root = ET.parse(svg_path).getroot()
    return float(root.get("width")), float(root.get("height"))


def _convert_svg(svg_path, out_path, ext):
    """R2DT only emits SVG; convert to whatever format --plot asked for via
    cairosvg. .png (the previous matplotlib-based renderer's only output
    format) and .pdf are supported; anything else falls back to .svg content
    copied as-is under the requested name, since silently writing nothing
    would be worse than an oddly-named SVG."""
    convert = _SVG_CONVERTERS.get(ext)
    if convert is None:
        logger.warning(f"--plot: {ext} not supported for R2DT output (only .svg, "
                       f".png, .pdf); writing raw SVG to {out_path} instead")
        shutil.copy(svg_path, out_path)
        return
    kwargs = {}
    if ext == ".png":
        width, height = _svg_intrinsic_size(svg_path)
        kwargs = {"scale": min(2.0, MAX_PNG_DIM / max(width, height, 1))}
    try:
        convert(url=svg_path, write_to=out_path, **kwargs)
    except Exception as e:
        logger.warning(f"SVG->{ext} conversion failed: {e}")


# --- main ---

def main():
    parser = argparse.ArgumentParser(
        description="assign Sprinzl coordinates to mt-tRNA sequences via structure-based cm selection.")
    parser.add_argument("--fasta", required=True,
                        help="input FASTA; headers: 'id|aa|anticodon|taxon', 'anticodon=XXX' tag, "
                             "or GtRNAdb-style 'tRNA-{AA}-{anticodon}' name (e.g. mt-tRNA-Ala-TGC-1-1)")
    parser.add_argument("--canonical-cm", required=True, nargs="+", metavar="CM_OR_DIR",
                        help="one or more canonical CM sources, tried in order per sequence: "
                             "a path to a single CM (e.g. TRNAinf-bact.cm, applies to every aa), "
                             "or a directory of {label}_{AA}.cm files (e.g. Metazoan_P.cm) to "
                             "select per-sequence by aa. the first source whose alignment anchors "
                             "the anticodon unambiguously is used; earlier sources take priority "
                             "(e.g. a bacterial CM first, then a metazoan per-AA directory, since "
                             "a CM built for the wrong clade can fail to thread a divergent loop)")
    parser.add_argument("--armless-cm-dir", required=True,
                        help="directory (searched recursively) for "
                             "armless_trn{AA}_wo_{d,t,d_and_t}.cm files")
    parser.add_argument("--out", default="sprinzl_mapping.tsv",
                        help="output TSV path (default: sprinzl_mapping.tsv)")
    parser.add_argument("--plot", default=None, metavar="SVG_OR_PNG_OR_PDF",
                        help="path for R2DT-rendered 2D diagram, one panel per "
                             "sequence arranged in a grid; format is chosen by "
                             "extension (.svg, .png, .pdf) (omit to skip plotting)")
    parser.add_argument("--ncols", type=int, default=6, help="plot grid columns")
    parser.add_argument("--r2dt-image", default=R2DT_DEFAULT_IMAGE, metavar="PATH",
                        help=f"R2DT Singularity image (default: {R2DT_DEFAULT_IMAGE})")
    parser.add_argument("-p", "--processes", type=int, default=4,
                        help="worker processes (default: 4)")
    parser.add_argument("--debug", action="store_true",
                        help="log alignment, arm-loss diagnosis, and CM routing for every sequence")
    args = parser.parse_args()

    if args.debug:
        _configure_logging("DEBUG")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    records = [(str(r.id) + (" " + r.description.split(None, 1)[1]
                              if " " in r.description else ""),
                str(r.seq))
               for r in SeqIO.parse(args.fasta, "fasta")]
    armless_cm_index = index_armless_cms(args.armless_cm_dir)

    canonical_cm_tiers = []
    tier_descs = []
    for source in args.canonical_cm:
        if os.path.isdir(source):
            tier_index = index_canonical_cms(source)
            canonical_cm_tiers.append(tier_index)
            tier_descs.append(f"per-AA from {source} ({len(tier_index)} CMs)")
        else:
            canonical_cm_tiers.append(source)
            tier_descs.append(source)

    logger.info(f"{len(records)} sequences, canonical CM tiers (in priority order): "
                f"{tier_descs}, "
                f"{len(armless_cm_index)} armless CMs available for rerouting, "
                f"{args.processes} worker process(es)")

    tasks = [(header, seq, canonical_cm_tiers, armless_cm_index, args.debug)
             for header, seq in records]

    if args.processes > 1:
        with multiprocessing.Pool(args.processes) as pool:
            results = pool.map(process_one_record, tasks)
    else:
        # single-process path: --debug logging interleaves in real time without
        # multiprocess log-buffering surprises.
        results = [process_one_record(t) for t in tasks]

    all_rows = [row for r in results for row in r["rows"]]
    n_failed = sum(1 for r in results if not r["rows"])
    if n_failed:
        logger.warning(f"{n_failed}/{len(records)} sequences produced no output")

    pd.DataFrame(all_rows).to_csv(args.out, sep="\t", index=False)
    logger.info(f"table: {args.out}")

    if args.plot:
        make_plot(results, args.plot, r2dt_image=args.r2dt_image, ncols=args.ncols)
        logger.info(f"plot:  {args.plot}")

        # sequences RNAfold-patched for a CM threading failure: also plot the
        # CM-only structure (pre-patch) side by side, so the patch's effect is
        # visible rather than assumed.
        cm_only_results = [{**r, "ss": r["cm_only_ss"]} for r in results if r.get("cm_only_ss")]
        if cm_only_results:
            cm_only_path = _cm_only_plot_path(args.plot)
            make_plot(cm_only_results, cm_only_path, r2dt_image=args.r2dt_image, ncols=args.ncols)
            logger.info(f"plot (CM-only, pre-RNAfold-patch): {cm_only_path}")


if __name__ == "__main__":
    main()


# =============================================================================
# arm-loss call string glossary
# every processed sequence emits exactly one of these in the log.
#
# CANONICAL_NO_ARM_LOSS -- every stem-loop slot passes absent() (see classify_arm_loss).
#
# T_OR_VAR_ARM_MISSING_slots=[n,..] -- one or more slots fail absent() (0-indexed,
#   5'->3'). a middle slot = optional variable-arm stem, never load-bearing on its
#   own (no armless CM for it). the LAST slot = T-arm: genuine loss, or a threading
#   failure patched via RNAfold -- see select_cm_and_align step 3.
#
# UPSTREAM_ARM_MISSING_offset=n -- D-arm absent, detected via register shift
#   (anticodon landed n slots downstream of expected). trusted directly, no
#   threading-failure cross-check (see select_cm_and_align).
#
# UPSTREAM_ARM_MISSING_slot=n -- D-arm absent with NO register shift (its own
#   slot n=idx-1 fails absent() instead). seen with CMs modeling more than the
#   canonical D/C/T trio. DOES get the threading-failure cross-check, unlike the
#   offset-based call above.
#
# BOTH_ARMS_MISSING_slots=[n,..] -- doubly-armless: D-arm and T-arm slots both
#   fail absent(), offset==0 (no single-arm shift with both arms gone). reroutes
#   to armless_trn{AA}_wo_d_and_t.cm.
#
# UNANCHORED_fallback_structurally_absent=[n,..] -- anticodon not uniquely
#   anchored (ambiguous AT-rich triplet); no directional signal, less reliable.
#
# threading failure (separate log line, not a call string): "CM diagnosed X-arm
#   missing (...) but the span folds as a real hairpin ... patching via RNAfold".
#   patch aborts silently (DEBUG log) on a bracket conflict with existing structure.
# =============================================================================
