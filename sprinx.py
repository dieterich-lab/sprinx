#!/usr/bin/env python3
"""
sprinx.py -- Sprinzl-coordinate annotation for mt-tRNAs.

1. problem
   mt-tRNAs exist in four structural shapes: cloverleaf, D-armless, T-armless,
   doubly-armless (Ozerova et al. 2024, PMC11571959). Sprinzl labels must be
   assigned relative to the correct shape; wrong shape -> wrong labels.

2. why score-based CM selection fails
   E-values are calibrated per model (Nawrocki & Eddy 2013, PMC3810854);
   an armless CM with fewer columns produces better E-values for canonical
   sequences than the canonical CM does, regardless of biological fit.
   length-normalising (bits/column) doesn't help: armless CMs retain the
   highest-information columns (acceptor + anticodon stems), inflating
   per-column scores. Rfam avoids this with hand-set per-family GA cutoffs
   (Kalvari et al., PMC6754622); this script avoids it by never comparing
   scores across models at all.

3. pipeline
   a. align to a canonical CM with cmalign --notrunc --nonbanded -g. --canonical-cm
      accepts multiple sources tried in priority order (e.g. bacterial whole-family
      CM, then a metazoan per-AA directory); the first whose anticodon anchors
      unambiguously wins -- never by score/E-value, see (2) -- because a CM built
      for the wrong clade can fail to thread a divergent sequence at all. details
      in select_cm_and_align.
   b. anchor on the anticodon; a missing UPSTREAM arm (D-arm) shifts remaining
      structure into wrong model columns (register shift). missing DOWNSTREAM
      arm (T-arm) does not shift. measure offset = expected_anticodon_slot - observed.
   c. n_pairs==0 at a stem slot means zero alignment columns have BOTH pairing
      partners simultaneously non-gap. this is a structural impossibility (if no
      column can form a pair, no stem can exist there), not a threshold judgment.
      however, n_pairs==0 has two distinct causes that require different responses:
        (i)  genuine arm loss: the sequence simply has no arm. the element span
             across the alignment is mostly or entirely gap characters.
        (ii) CM threading failure: the arm exists but cmalign placed its sequence
             into unmodeled insert columns because the arm is too divergent from
             the CM consensus. the stem model columns are all gaps, but the span
             DOES contain nucleotides as insert characters.
      distinguishing (i) from (ii): count non-gap nucleotides across the full
      element span (stem + loop model columns + intervening insert characters).
      if the count is < n_stem_cols + MIN_HAIRPIN_LOOP (=3, steric minimum for
      RNA backbone; Woodson 2008, doi:10.1146/annurev.physchem.59.032607.093743),
      no hairpin can form physically: genuine arm loss (i). otherwise: threading
      failure (ii).
      hybrid Infernal + RNAfold design: for threading failures, Infernal's
      canonical CM is correct for all arms it DID thread properly; only the
      mis-threaded arm needs structural recovery. RNAfold MFE on the short arm
      span (typically 13-20 nt) is reliable at this length because competing folds
      are energetically negligible. the hybrid avoids two failure modes: (a)
      relying on Infernal alone would call threading failures as arm loss and
      misroute to an armless CM; (b) relying on RNAfold alone for full-sequence
      mt-tRNA folding is unreliable due to tertiary interactions and base
      modifications not captured by 2D MFE (Helm 2006, doi:10.1093/nar/gkl348).
   d. if genuinely absent, reroute to armless CM (Ozerova et al. 2024,
      PMC11571959). isoacceptors (Leu1/Leu2, Ser1/Ser2) disambiguated by
      anticodon, not filename suffix. for doubly-armless (D + T both missing),
      routes to the d_and_t CM.
   e. assign Sprinzl coordinates (Sprinzl et al. 1998, PMC147216).

4. implementation notes
   cmalign flags (required together, every call):
     --notrunc   : include all positions; without it, local mode silently drops
                   regions that fit poorly, causing false arm-loss calls.
     --nonbanded : exact CYK/Inside DP; HMM banding is ~10x faster but
                   introduces alignment errors on divergent mt-tRNA structures.
     -g          : glocal; prevents local begin/end states skipping arm regions.
     flag set follows QutRNA2 (github.com/dieterich-lab/QutRNA2,
     biorxiv:2025.10.20.683443).
   header format (pipe-delimited):
     field 1: seq id | field 2: three-letter aa (e.g. Ala, Leu1)
     field 3: anticodon (3nt, RNA or DNA) | field 4: taxon
     fallback: 'anticodon=XXX' tag anywhere in the header.
     field 3 is the primary key for CM selection; field 2 only identifies aa.
   armless CM filenames: armless_trn{AA}_wo_{arm}.cm where arm is d, t, or d_and_t
   for doubly-armless (Ozerova et al. 2024, PMC11571959). armless CM rerouting
   is unaffected by which canonical CM tier won above -- it only triggers once
   a genuine arm-loss diagnosis is made from whichever tier's alignment was used.
   each --canonical-cm source is a directory of {label}_{AA}.cm files (e.g.
   Metazoan_P.cm; label/clade is ignored, selection is by AA only, per-sequence,
   same as armless CM selection) or a single CM file (applies to every aa,
   e.g. a whole-family CM like TRNAinf-bact.cm).

5. output
   sprinzl_mapping.tsv: seq_id, seq_index, nucleotide, sprinzl_position, region,
   cm_used, rerouted, arm_loss_call. optional --plot for 2D cloverleaf PNGs.

usage:
  python sprinx.py --fasta seqs.fa --canonical-cm TRNAinf-euk.cm \\
      --armless-cm-dir cm_models/ --out results/sprinzl_mapping.tsv
  python sprinx.py --fasta seqs.fa --canonical-cm TRNAinf-euk.cm \\
      --armless-cm-dir cm_models/ --plot --processes 8 --debug
  python sprinx.py --fasta seqs.fa --canonical-cm cm_models_by_clade/ \\
      --armless-cm-dir cm_models/ --processes 8
  # multiple canonical CM tiers, tried in order (bacterial whole-family CM first,
  # then a metazoan per-AA directory)
  python sprinx.py --fasta seqs.fa \\
      --canonical-cm TRNAinf-bact.cm cm_models_by_clade/ \\
      --armless-cm-dir cm_models/ --processes 8
"""

import argparse
import multiprocessing
import os
import re
import subprocess
import sys
import tempfile
import warnings
from collections import defaultdict

import matplotlib
import numpy as np
import pandas as pd
import RNA
from forgi.graph.bulge_graph import BulgeGraph
from Bio.Data.IUPACData import protein_letters_3to1
from Bio import SeqIO
from loguru import logger
from scipy.stats import binomtest

import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
matplotlib.use("Agg")  # must precede any pyplot state; no display available

warnings.filterwarnings("ignore")

def _configure_logging(level):
    """(re)point loguru at stderr with a bare message format. called at import
    (INFO) and again wherever --debug is honoured (main(), each worker process
    since multiprocessing forks/spawns fresh interpreters)."""
    logger.remove()
    logger.add(sys.stderr, format="<level>{message}</level>", level=level)


_configure_logging("INFO")


# --- constants: tRNA topology facts + Sprinzl coordinate system (PMC147216) ---

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
for _p in range(49, 54):  SPRINZL_REGION[str(_p)] = "T_stem_5"
for _p in range(54, 61):  SPRINZL_REGION[str(_p)] = "T_loop"
for _p in range(61, 66):  SPRINZL_REGION[str(_p)] = "T_stem_3"
for _p in range(66, 73):  SPRINZL_REGION[str(_p)] = "acceptor_3"
for _p in range(73, 77):  SPRINZL_REGION[str(_p)] = "discriminator_CCA"

NT_COLOUR = {"A": "#c0392b", "U": "#2471a3", "G": "#d35400", "C": "#1e8449"}

# armless CM filename regex; naming follows Ozerova et al. 2024 (PMC11571959).
ARMLESS_CM_RE = re.compile(r"armless_trn(\w+)_wo_(d_and_t|d|t)\.cm$")


# --- generic helpers ---

def run(cmd):
    """run a subprocess; log the exact command at debug level for manual reproduction."""
    logger.debug(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    return result.stdout, result.stderr, result.returncode


def pair_table(ss):
    """dot-bracket -> {i: j, j: i}, 0-indexed. uses RNA.ptable() which raises
    on malformed input instead of silently mishandling an unpartnered bracket;
    callers that may produce orphan brackets should run drop_orphan_brackets first."""
    pt = RNA.ptable(ss)
    if pt is None:
        raise ValueError(f"ViennaRNA rejected this as a valid structure: {ss!r}")
    return {i - 1: pt[i] - 1 for i in range(1, len(pt)) if pt[i] > 0}


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


def header_to_anticodon(header):
    """extract anticodon from 'id|aa|anticodon|taxon' (field 3) or 'anticodon=XXX'
    tag anywhere in the header. returns 3-nt RNA string or None; returns None
    rather than guessing on format mismatch -- a wrong anticodon propagates
    through the entire Sprinzl assignment."""
    fields = header.split("|")
    if len(fields) >= 3 and re.fullmatch(r"[ACGUTacgut]{3}", fields[2]):
        return fields[2].upper().replace("T", "U")
    m = re.search(r"anticodon=([ACGUTacgut]{3})", header)
    return m.group(1).upper().replace("T", "U") if m else None


def aa_field_to_cm_code(aa_field, cm_index_keys):
    """header aa field (e.g. 'Ala', 'Leu1') -> one-letter CM code (e.g. 'A', 'L1').
    derivation: strip digit suffix, protein_letters_3to1 (IUPAC 1984), reattach
    suffix, check against cm_index_keys. returns None if code absent from index.
    isoacceptor disambiguation uses the anticodon (resolve_armless_cm), not this."""
    if not aa_field:
        return None
    m = re.fullmatch(r"([A-Za-z]+)(\d*)", aa_field.strip())
    if not m:
        return None
    one = protein_letters_3to1.get(m.group(1).capitalize())
    if one is None:
        return None
    code = one + m.group(2)
    return code if code in {aa for aa, _ in cm_index_keys} else None


def header_to_aa(header):
    """return the raw aa field (second pipe-delimited field) or None."""
    fields = header.split("|")
    return fields[1].strip() if len(fields) >= 2 and fields[1].strip() else None


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


def _resolve_canonical_for_tier(header, tier):
    """resolve one canonical-CM tier to a concrete .cm path for this header, or
    None if the tier doesn't apply (a per-AA dict with no entry for this aa).
    a plain path string tier applies unconditionally (e.g. a whole-family CM
    like TRNAinf-bact.cm, which models every amino acid with one CM)."""
    if isinstance(tier, dict):
        aa_code = aa_field_to_cm_code(header_to_aa(header), {(aa, None) for aa in tier})
        return tier.get(aa_code)
    return tier


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
    numbering. strips both gap symbols: '-' (model deletion) and '.' (multi-seq
    insert gap). converts WUSS to dot-bracket; repairs orphan brackets from
    gap-stripping via drop_orphan_brackets."""
    aligned_seq, ss_cons = alignment["aligned_seq"], alignment["ss_cons"]
    db = RNA.db_from_WUSS(ss_cons)
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

def get_stem_loop_elements(ss):
    """decompose dot-bracket (or WUSS) into ordered inner stem-loop dicts,
    merging across internal-loop bulges via forgi BulgeGraph. excludes the outer
    acceptor stem (no hairpin loop). each dict: {'stem_cols', 'loop_cols', 'span'},
    all 0-indexed into the alignment (or sequence if ss is ungapped).
    RNA.db_from_WUSS is idempotent on plain dot-bracket, so this accepts either."""
    db = RNA.db_from_WUSS(ss)
    bg = BulgeGraph.from_dotbracket(db)
    stem_elems = sorted([e for e in bg.defines if e.startswith("s")],
                        key=lambda e: bg.defines[e][0])

    def interior_neighbors(elem):
        # stem elements on the other side of a single interior loop (bulge merging)
        return {nb2 for nb in bg.edges[elem] if nb.startswith("i")
                for nb2 in bg.edges[nb] if nb2.startswith("s") and nb2 != elem}

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
        used |= group

        stem_cols = sorted(set(
            col for g in group
            for a, b, c, d in [bg.defines[g]]
            for col in list(range(a - 1, b)) + list(range(c - 1, d))
        ))
        loop_cols = []
        for g in group:
            for nb in bg.edges[g]:
                if nb.startswith("h"):
                    a, b = bg.defines[nb][:2]
                    loop_cols = list(range(a - 1, b))

        groups.append({
            "stem_cols": stem_cols,
            "loop_cols": loop_cols,
            "span": (min(stem_cols), max(stem_cols) + 1),
        })

    groups.sort(key=lambda g: g["span"][0])
    return [g for g in groups if g["loop_cols"]]  # drop acceptor stem (no hairpin)


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
    """WC/wobble pairing check for one stem element.
    n_pairs: columns where both pairing partners are simultaneously non-gap.
    n_compatible: of those, pairs that are WC or G-U wobble.
    p_value: binomial test against null rate len(WC_PAIRS)/16 = 6/16 (the fraction
      of the 16 possible dinucleotides that form a WC/wobble pair; this is a
      combinatorial fact, not a fitted parameter).
    n_pairs==0 is a statement of structural impossibility, not a threshold call:
      if every alignment column that the CM designated as a stem pair has at least
      one of its two partners as a gap character, then no base pair can form at any
      position in the stem regardless of the nucleotides present. the arm cannot
      exist. this is categorically different from n_pairs > 0 with low n_compatible,
      which would suggest a poorly-paired but structurally present stem.
    p_value is not thresholded: short stems (3-5 bp, common for D-arm) frequently
      lack statistical power even when the stem is real. callers read
      per_stem_complementarity directly rather than relying on a binary verdict.
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
# triloops are the smallest observed RNA hairpins. Woodson 2008, Annu Rev Phys
# Chem 59:617-639. doi:10.1146/annurev.physchem.59.032607.093743
MIN_HAIRPIN_LOOP = 3

# soft threshold used by classify_arm_loss's absent(): 1-2 coincidental base
# pairs are too few to nucleate a stable helix, so n_pairs<3 is weak evidence
# a stem is real (3 is empirically the smallest count with no false positives
# on the canonical-36 test set). UNLIKE MIN_HAIRPIN_LOOP (geometric certainty),
# this is a judgment call, not a certainty -- it only flags *candidates* for
# arm loss; every candidate still has to pass the hard arm_span_has_enough_sequence
# check before any reroute happens.
MIN_STEM_PAIRS = 3


def arm_span_has_enough_sequence(aligned_seq, elem):
    """first-stage (fast, hard) filter after a stem slot is flagged absent: does
    the span even contain enough nucleotides to physically form a hairpin
    (n_stem_cols + MIN_HAIRPIN_LOOP, the steric minimum; Woodson 2008,
    doi:10.1146/annurev.physchem.59.032607.093743)? below that: definite genuine
    loss. at or above: NOT proof the arm is present, just not ruled out by
    volume alone -- a CM with wide insert-state capacity can pass this on
    unrelated leftover sequence alone (see arm_is_threading_failure, the
    required second-stage check that actually folds the span).
    returns False for definite genuine loss, True to proceed to that check."""
    start, end = elem["span"]
    n_nts = sum(1 for c in aligned_seq[start:end] if c not in "-.")
    return n_nts >= len(elem["stem_cols"]) + MIN_HAIRPIN_LOOP


def _arm_insert_subseq_and_fold(aligned_seq, final_seq, elem):
    """extract ONLY the mis-threaded insert-state (lowercase) nucleotides within
    elem's span and fold them with RNAfold MFE. narrow on purpose -- used for
    PATCHING (patch_threading_failure_arm), so already-correct matched columns
    are never touched. NOT for detecting whether a hairpin exists at all; a
    real arm isn't always in insert columns -- see _arm_full_span_fold for that.
    returns (ungapped_positions, arm_ss), or (None, None) if there's too
    little sequence to fold (< MIN_HAIRPIN_LOOP + 2 nt)."""
    span_start, span_end = elem["span"]
    ungapped_positions, ungapped_idx = [], 0
    for gapped_idx, c in enumerate(aligned_seq):
        if c not in "-.":
            if span_start <= gapped_idx < span_end and c.islower():
                ungapped_positions.append(ungapped_idx)
            ungapped_idx += 1

    if len(ungapped_positions) < MIN_HAIRPIN_LOOP + 2:
        return None, None

    arm_subseq = "".join(final_seq[p] for p in ungapped_positions)
    arm_ss, _ = RNA.fold_compound(arm_subseq).mfe()
    return ungapped_positions, arm_ss


def _arm_full_span_fold(aligned_seq, final_seq, elem):
    """fold the FULL non-gap span (matched columns + inserts together), unlike
    _arm_insert_subseq_and_fold's insert-only extraction -- used by
    arm_is_threading_failure to detect a hairpin anywhere in the span.
    returns arm_ss, or None if there's too little sequence to fold
    (< MIN_HAIRPIN_LOOP + 2 nt)."""
    span_start, span_end = elem["span"]
    ungapped_positions, ungapped_idx = [], 0
    for gapped_idx, c in enumerate(aligned_seq):
        if c not in "-.":
            if span_start <= gapped_idx < span_end:
                ungapped_positions.append(ungapped_idx)
            ungapped_idx += 1

    if len(ungapped_positions) < MIN_HAIRPIN_LOOP + 2:
        return None

    arm_subseq = "".join(final_seq[p] for p in ungapped_positions)
    arm_ss, _ = RNA.fold_compound(arm_subseq).mfe()
    return arm_ss


def arm_is_threading_failure(aligned_seq, final_seq, elem):
    """second-stage check, run only after arm_span_has_enough_sequence passes:
    does the span actually fold as a hairpin? that raw-count check can be
    fooled by a CM with wide insert-state capacity (e.g. a whole-family
    bacterial CM) -- a genuinely absent arm's span can pass on volume alone,
    filled with unrelated leftover sequence (e.g. 3' trailer) that folds as
    nothing. folding the whole span here (both matched and insert columns,
    since a real arm isn't always mis-threaded into inserts -- see
    _arm_full_span_fold) is a much stronger positive signal.
    True: RNAfold found a hairpin -- real, recoverable arm.
    False: nothing folds -- genuine loss despite passing the count check."""
    arm_ss = _arm_full_span_fold(aligned_seq, final_seq, elem)
    return arm_ss is not None and "(" in arm_ss


def patch_threading_failure_arm(aligned_seq, final_seq, final_ss, elem):
    """recover arm structure for a CM threading failure: arm present but mis-threaded
    into insert columns. called only once arm_is_threading_failure has confirmed
    a real hairpin exists there. the rest of final_ss (from the canonical CM) is
    already correct, so we splice in ONLY the mis-threaded span's own RNAfold
    fold rather than refolding the whole molecule -- full-sequence RNAfold on a
    mt-tRNA is unreliable (tertiary contacts, base modifications; Helm 2006,
    doi:10.1093/nar/gkl348), but a short isolated span (13-20nt) is fine.
    safety: every bracket in the fold must target a '.' in final_ss, or the
    patch is aborted -- placing an open whose matching close is blocked leaves
    a dangling bracket RNA.ptable() would reject.
    returns patched final_ss, or the original if: no stem found, a bracket
    conflicts with existing structure, or the patch is unbalanced (safety net)."""
    ungapped_positions, arm_ss = _arm_insert_subseq_and_fold(aligned_seq, final_seq, elem)

    if ungapped_positions is None or "(" not in arm_ss:
        return final_ss

    # pre-validate: every bracket in arm_ss must target a '.' in final_ss.
    # a rejected close leaves its matching open dangling -> unbalanced.
    ss_list = list(final_ss)
    for i, c in enumerate(arm_ss):
        if c in "()" and i < len(ungapped_positions):
            if ss_list[ungapped_positions[i]] != ".":
                logger.debug(
                    f"arm patch aborted: arm_ss[{i}]={c!r} targets "
                    f"final_ss[{ungapped_positions[i]}]="
                    f"{ss_list[ungapped_positions[i]]!r}; span overlaps existing structure"
                )
                return final_ss

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


def resolve_armless_cm(header, seq, aa_code, missing_arm, anticodon, armless_cm_index):
    """pick the correct armless CM for this sequence and arm type. if multiple
    isoacceptor CMs share the same base aa code (Leu1/Leu2, Ser1/Ser2), align to
    each and return the one where the header anticodon lands in the anticodon loop.
    the anticodon is the discriminating fact, not filename suffix. falls back to
    first candidate on disambiguation failure; returns None if no CM exists."""
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
    if anticodon:
        for path in candidates:
            aln = cmalign_one(header, seq, path)
            if aln is None:
                continue
            elements = get_stem_loop_elements(aln["ss_cons"])
            if find_anticodon_stem_index(aln["aligned_seq"], elements, anticodon)[0] is not None:
                logger.debug(f"{header}: isoacceptor {anticodon} -> {os.path.basename(path)}")
                return path
    logger.warning(f"{header}: anticodon disambiguation failed for aa={aa_code!r} arm={missing_arm!r}; "
                   f"using {os.path.basename(candidates[0])}")
    return candidates[0]


def _routing_result(final_alignment, cm_used, diagnosis, rerouted=False, threading_failure_elem=None):
    """assemble the dict select_cm_and_align returns at each of its exit points,
    so the shape is defined once instead of copy-pasted per branch."""
    return {"final_alignment": final_alignment, "cm_used": cm_used, "diagnosis": diagnosis,
            "rerouted": rerouted, "threading_failure_elem": threading_failure_elem}


def select_cm_and_align(header, seq, canonical_cm_tiers, armless_cm_index):
    """top-level CM selection for one sequence. routing order:
      1 try each canonical CM tier in order (e.g. bacterial whole-family CM,
        then a metazoan per-AA directory); align + diagnose against each and
        take the first whose anticodon anchors unambiguously. a tier that
        doesn't apply to this aa (per-AA dict with no matching entry) or
        whose cmalign fails outright is skipped. if no tier anchors cleanly,
        fall back to the first tier that produced any alignment at all --
        same graceful degradation as the old single-canonical-CM behaviour.
        this exists because a CM built for the wrong clade can lack the
        capacity to model a divergent loop (e.g. an unusually long variable
        loop), causing cmalign to thread the overflow into an adjacent arm's
        insert states and break anticodon anchoring entirely -- a different
        CM (e.g. bacterial, given mitochondria's endosymbiotic origin) may
        thread the same sequence correctly. selection is never by alignment
        score/E-value across tiers, only by whether the anchor is clean --
        see module docstring section 2 for why cross-model score comparison
        is invalid.
      2 if no arm missing (or only variable arm): return canonical alignment.
      3 if D-arm (no register shift) or T-arm flagged absent (classify_arm_loss's
        absent(), a soft MIN_STEM_PAIRS-based signal): two checks, both required
        to call it a threading failure rather than genuine loss --
        arm_span_has_enough_sequence (fast, hard: enough nt to physically form
        a hairpin?) then arm_is_threading_failure (RNAfold: does it actually
        fold as one?). failing either -> genuine loss, proceed to rerouting.
        passing both -> keep canonical alignment, patch via RNAfold instead.
        D-arm loss found via register shift (offset>0) is trusted directly and
        skips both checks -- see the comment below.
      4 reroute: resolve_armless_cm with anticodon-based isoacceptor disambiguation.
      5 no matching armless CM: warn and return canonical alignment.
    canonical_cm_tiers is normally a list of tiers (each a .cm path applying to
    every sequence, or a {aa_code: path} dict from index_canonical_cms), tried
    in order; a bare path string or dict is also accepted and wrapped as a
    single-tier list, for callers (and tests) that only have one CM.
    returns dict: final_alignment, cm_used, diagnosis, rerouted, threading_failure_elem."""
    if isinstance(canonical_cm_tiers, (str, dict)):
        canonical_cm_tiers = [canonical_cm_tiers]

    canonical_alignment = canonical_cm = diagnosis = None
    for tier in canonical_cm_tiers:
        path = _resolve_canonical_for_tier(header, tier)
        if path is None:
            continue
        aln = cmalign_one(header, seq, path)
        if aln is None:
            continue
        diag = classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
        if canonical_alignment is None:               # keep first usable tier as fallback
            canonical_alignment, canonical_cm, diagnosis = aln, path, diag
        if diag["anticodon_stem_index"] is not None:   # clean anchor: stop searching
            canonical_alignment, canonical_cm, diagnosis = aln, path, diag
            break

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
                    f"{header}: CM diagnosed {missing_arm}-arm missing ({diagnosis['call']}) "
                    f"but the span folds as a real hairpin "
                    f"(CM threading failure, not genuine arm loss) -- patching via RNAfold"
                )
                return _routing_result(canonical_alignment, canonical_cm, diagnosis,
                                        threading_failure_elem=elem)

    # genuine arm loss: reroute
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
    """acceptor stem (first '(' run, capped at 7bp) + inner stems. does not
    label D/C/T yet -- that needs the anticodon (see locate_anticodon_stem)."""
    pairs = pair_table(ss)
    n = len(ss)
    start = next((i for i in range(n) if ss[i] == "("), None)
    if start is None:
        raise ValueError(f"no acceptor stem (no '(' at all) in: {ss!r}")
    acceptor = []
    i = start
    while i < n and ss[i] == "(" and len(acceptor) < 7:
        acceptor.append(i)
        i += 1
    acceptor_pairs = sorted(pairs[p] for p in acceptor)
    inner_start, inner_end = acceptor[-1] + 1, acceptor_pairs[0]

    inner_stems, i = [], inner_start
    while i < inner_end:
        if ss[i] == "(":
            stem5 = []
            while i < inner_end and ss[i] == "(":
                stem5.append(i)
                i += 1
            inner_stems.append(stem5)
        else:
            i += 1

    return {
        "acceptor_5": acceptor, "acceptor_3": acceptor_pairs,
        "inner_stems": inner_stems, "inner_start": inner_start, "inner_end": inner_end,
        "trailer": [p for p in range(acceptor_pairs[-1] + 1, n) if ss[p] == "."],
    }


def locate_anticodon_stem(topo, ss, seq, anticodon):
    """identify C-stem (anticodon arm) by anticodon content, D-stem = sibling
    before C that does not enclose it, T-stem = sibling after it.
    innermost-stem-first search order prevents a D-armless pseudostem (which
    structurally encloses C) from being matched before the real C-stem.
    'does not enclose' check: a pseudostem that opens before C and closes after
    C must be excluded as D-arm candidate -- it IS the enclosing pseudostem.
    see TestSprinzlAssignment::test_d_armless_replacement_loop_gets_d_arm_labels."""
    pairs = pair_table(ss)
    inner_stems = topo["inner_stems"]

    def full_loop(stem5):
        return list(range(stem5[-1] + 1, pairs[stem5[-1]])) if stem5 else []

    def direct_loop(idx, stem5):
        # loop positions exclusive of any nested stem's span
        lp = set(full_loop(stem5))
        for j, other5 in enumerate(inner_stems):
            if j != idx and stem5[-1] < other5[0] < pairs[stem5[-1]]:
                lp -= set(range(other5[0], pairs[other5[-1]] + 1))
        return sorted(lp)

    def unpaired(a, b):
        return [p for p in range(a, b) if ss[p] == "."]

    ac = (anticodon or "").upper().replace("T", "U")
    search_order = sorted(range(len(inner_stems)),
                          key=lambda i: pairs[inner_stems[i][-1]] - inner_stems[i][0])

    c_stem, c_idx = [], None
    for idx in search_order:
        stem5 = inner_stems[idx]
        if ac and ac in "".join(seq[p] for p in direct_loop(idx, stem5)):
            c_stem, c_idx = stem5, idx
            break
    if c_stem is None and inner_stems:
        c_idx = search_order[0]
        c_stem = inner_stems[c_idx]

    c_close = pairs[c_stem[-1]] if c_stem else None
    d_stem, t_stem = [], []
    if c_stem:
        before = [s5 for s5 in inner_stems
                  if s5 != c_stem and s5[0] < c_stem[0] and pairs[s5[-1]] < c_close]
        after = [s5 for s5 in inner_stems if s5 != c_stem and s5[0] > c_close]
        d_stem = max(before, key=len) if before else []
        # t-arm is the last (highest-position) stem after c-close.
        # min(after) breaks for class-ii tRNAs (ser, leu) and some mt-tRNAs
        # with a variable arm stem: it picks the variable arm as t-arm instead.
        t_stem = max(after, key=lambda s5: s5[0]) if after else []

    d_close = pairs[d_stem[-1]] if d_stem else None
    t_open = t_stem[0] if t_stem else None

    if c_close is not None and t_open is not None:
        # all positions between c-stem close and t-arm open, not filtered to
        # unpaired: class-ii tRNAs have a variable arm stem in this region
        # whose paired positions must also receive sprinzl 44-48 labels.
        var_loop = list(range(c_close + 1, t_open))
    elif c_close is not None:
        var_loop = list(range(c_close + 1, topo["inner_end"]))
    else:
        var_loop = []

    linker_5_end = d_stem[0] if d_stem else (c_stem[0] if c_stem else topo["inner_end"])

    return {
        "d_stem5": d_stem, "d_stem3": sorted(pairs[p] for p in d_stem) if d_stem else [],
        "d_loop": full_loop(d_stem),
        "c_stem5": c_stem, "c_stem3": sorted(pairs[p] for p in c_stem) if c_stem else [],
        "c_loop": full_loop(c_stem),
        "t_stem5": t_stem, "t_stem3": sorted(pairs[p] for p in t_stem) if t_stem else [],
        "t_loop": full_loop(t_stem),
        "var_loop": var_loop,
        "linker_5": unpaired(topo["inner_start"], linker_5_end),
        "linker_dc": unpaired(d_close + 1, c_stem[0]) if (d_close is not None and c_stem) else [],
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


def sprinzl_map(ss, seq, anticodon):
    """assign a Sprinzl label to every nucleotide index; returns {seq_index: label}.
    D-armless tRNAs: replacement loop (all of linker_5) is mapped onto D-arm Sprinzl
    positions 8-26 by structural analogy, following Ozerova et al. 2024 (PMC11571959).
    missing T-arm produces no labels for its region."""
    topo = parse_topology(ss)
    arms = locate_anticodon_stem(topo, ss, seq, anticodon)
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
    assign_slots(labels, arms["c_loop"],     [str(i) for i in range(32, 39)])
    assign_slots(labels, arms["c_stem3"],    [str(i) for i in range(39, 44)])
    assign_slots(labels, arms["var_loop"],   ["44", "45", "46", "47", "48"])
    assign_slots(labels, arms["t_stem5"],    [str(i) for i in range(49, 54)])
    assign_slots(labels, arms["t_loop"],     [str(i) for i in range(54, 61)])
    assign_slots(labels, arms["t_stem3"],    [str(i) for i in range(61, 66)])
    return labels


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
    if routing.get("threading_failure_elem"):
        final_ss = patch_threading_failure_arm(
            alignment["aligned_seq"], final_seq, final_ss,
            routing["threading_failure_elem"]
        )

    if len(final_seq) != len(seq):
        logger.warning(f"{header}: ungapped length {len(final_seq)} != input {len(seq)}, skipped")
        return {"header": header, "rows": [], "summary": "LENGTH_MISMATCH"}

    anticodon = header_to_anticodon(header)
    if anticodon is None:
        logger.warning(f"{header}: no anticodon in header; C-stem location unreliable")

    sprinzl = sprinzl_map(final_ss, final_seq, anticodon)

    diagnosis = routing["diagnosis"] or {}
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
        region_key = re.match(r"\d+", label).group() if label else ""
        rows.append({
            "seq_id": header, "seq_index": i, "nucleotide": base,
            "sprinzl_position": label, "region": SPRINZL_REGION.get(region_key, ""),
            "cm_used": cm_name, "rerouted": routing["rerouted"],
            "arm_loss_call": diagnosis.get("call"),
        })

    logger.info(f"{header}: {summary}  [{diagnosis.get('call')}]")
    return {"header": header, "rows": rows, "summary": summary,
            "seq": final_seq, "ss": final_ss, "sprinzl": sprinzl}


# --- plotting (optional, --plot) ---

def naview_coords(ss):
    """2D coordinates from ViennaRNA NAVIEW layout algorithm."""
    RNA.cvar.rna_plot_type = 1
    obj = RNA.get_xy_coordinates(ss)
    return np.array([[obj.get(i).X, obj.get(i).Y] for i in range(len(ss))])


def orient_acceptor_north(coords, acceptor_idx, anticodon_loop_idx):
    """rotate layout so acceptor stem points up, anticodon loop points down."""
    if not acceptor_idx or not anticodon_loop_idx:
        return coords
    vec = coords[acceptor_idx].mean(0) - coords[anticodon_loop_idx].mean(0)
    angle = np.pi / 2 - np.arctan2(vec[1], vec[0])
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    centre = coords.mean(0)
    return (coords - centre) @ rot.T + centre


def draw_trna(ax, seq, ss, sprinzl, title, label_step=5, fontsize=8):
    """backbone, base-pair lines, nucleotide letters, and Sprinzl labels on ax."""
    pairs = pair_table(ss)
    topo = parse_topology(ss)
    coords = naview_coords(ss)
    anticodon_loop = []
    for stem5 in topo["inner_stems"]:
        anticodon_loop = list(range(stem5[-1] + 1, pairs[stem5[-1]]))
    coords = orient_acceptor_north(coords, topo["acceptor_5"] + topo["acceptor_3"], anticodon_loop)

    for i in range(len(seq) - 1):
        ax.plot(*coords[[i, i + 1]].T, "-", color="#c0c0c0", lw=0.8, zorder=1)
    for i, j in pairs.items():
        if i < j:
            ax.plot(*coords[[i, j]].T, "-", color="#909090", lw=0.7, alpha=0.7, zorder=1)
    for i, base in enumerate(seq):
        ax.text(*coords[i], base, ha="center", va="center", fontsize=fontsize,
                fontweight="bold", color=NT_COLOUR.get(base, "#555555"),
                fontfamily="monospace", zorder=3)
        label = sprinzl.get(i, "")
        # label every Nth integer position; always show lettered inserts (17a, 60A, ...)
        show = label and (
            not label[:-1].isdigit()
            or int(re.match(r"\d+", label).group()) % label_step == 0
            or label == "1"
        )
        if show:
            ax.text(coords[i, 0], coords[i, 1] - 2.3, label, ha="center", va="top",
                    fontsize=fontsize * 0.6, color="#333333", zorder=4)

    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=fontsize + 2, pad=6)
    margin = 5
    ax.set_xlim(coords[:, 0].min() - margin, coords[:, 0].max() + margin)
    ax.set_ylim(coords[:, 1].min() - margin, coords[:, 1].max() + margin)


def make_plot(results, out_path, ncols=4, label_step=5):
    """grid of one cloverleaf per successfully-processed record; failed records skipped.
    plotted in header order: species (taxon field) first, then tRNA (aa field), so
    isoacceptors of the same species group together and species cluster in the grid."""
    plotted = [r for r in results if r["rows"]]
    if not plotted:
        logger.warning("nothing to plot -- every record failed upstream")
        return
    plotted.sort(key=lambda r: (header_to_taxon(r["header"]) or "",
                                header_to_aa(r["header"]) or "", r["header"]))
    ncols = min(ncols, len(plotted))
    nrows = -(-len(plotted) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 5, nrows * 6))
    axes = np.atleast_1d(axes).flatten()
    for ax, r in zip(axes, plotted):
        draw_trna(ax, r["seq"], r["ss"], r["sprinzl"],
                  title=f"{r['header']}\n{r['summary']}", label_step=label_step)
    for ax in axes[len(plotted):]:
        ax.axis("off")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# --- main ---

def main():
    parser = argparse.ArgumentParser(
        description="assign Sprinzl coordinates to mt-tRNA sequences via structure-based cm selection.")
    parser.add_argument("--fasta", required=True,
                        help="input FASTA; headers: 'id|aa|anticodon|taxon' or 'anticodon=XXX' tag")
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
    parser.add_argument("--plot", default=None, metavar="PNG",
                        help="path for PNG cloverleaf grid (omit to skip plotting)")
    parser.add_argument("--ncols", type=int, default=4, help="plot grid columns")
    parser.add_argument("--label-step", type=int, default=5,
                        help="plot: label every Nth Sprinzl position")
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
        make_plot(results, args.plot, ncols=args.ncols, label_step=args.label_step)
        logger.info(f"plot:  {args.plot}")


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
