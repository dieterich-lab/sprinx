"""
sprinx.mito: canonical-CM tiering, arm-loss diagnosis, and armless-CM
rerouting for mt-tRNAs.

Everything here exists because mt-tRNAs can genuinely lose an arm (D-arm,
T-arm, or both - Ozerova et al. 2024); cytosolic/nuclear tRNAs (sprinx.cyto)
don't have this problem, so none of this module applies to them. Structural
parsing and Sprinzl-label assignment (forgi topology, sprinzl_map) are
generic and live in sprinx.common instead.

1. why score-based CM selection fails
   E-values are calibrated per model (Infernal User Guide); an armless CM
   with fewer columns produces better E-values for canonical sequences than
   the canonical CM does, regardless of biological fit. length-normalising
   (bits/column) doesn't help: armless CMs retain the highest-information
   columns (acceptor + anticodon stems), inflating per-column scores. Rfam
   avoids this with hand-set per-family GA cutoffs; this module avoids it by
   never comparing scores across models of different structure at all - only
   across canonical tiers, scored by total base-pairing evidence (see
   select_cm_and_align).

2. pipeline
   a. align to a canonical CM with cmalign --notrunc --nonbanded -g. --canonical-cm
      accepts multiple sources tried in priority order (e.g. bacterial whole-family
      CM, then a metazoan per-AA directory); among tiers that anchor the
      anticodon, the one with the highest total base-pairing evidence wins.
      details in select_cm_and_align.
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
   e. assign Sprinzl coordinates (sprinx.common.sprinzl_map).

3. armless CM filenames: armless_trn{AA}_wo_{arm}.cm where arm is d, t, or
   d_and_t for doubly-armless (Ozerova et al. 2024). armless CM rerouting
   is unaffected by which canonical CM tier won above; it only triggers once
   an arm-loss diagnosis is made from whichever tier's alignment was used.
   each --canonical-cm source is a directory of {label}_{AA}.cm files (e.g.
   Metazoan_P.cm; label/clade is ignored, selection is by AA only, per-sequence,
   same as armless CM selection) or a single CM file (applies to every aa,
   e.g. a whole-family CM like TRNAinf-bact.cm).

4. output
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

import RNA
from loguru import logger

from sprinx.common import (
    SPRINZL_REGION,
    _configure_logging,
    _forgi_stem_groups,
    _pick_by_anticodon_anchor,
    _scan_cm_files,
    aa_field_to_cm_code,
    cmalign_one,
    finalize_structure,
    find_anticodon_stem_index,
    get_stem_loop_elements,
    header_to_aa,
    header_to_anticodon,
    package_data_path,
    slide_stems_to_improve_pairing,
    sprinzl_map,
    stem_complementarity,
)


def default_canonical_cm_sources():
    """bundled default --canonical-cm sources: bacterial whole-family CM
    first, then a per-AA metazoan directory. covers metazoan mitochondrial
    tRNAs; a different clade needs its own CMs supplied via --canonical-cm."""
    canonical_dir = package_data_path("mito_cm", "canonical")
    return [
        os.path.join(canonical_dir, "TRNAinf-bact.cm"),
        os.path.join(canonical_dir, "mitofinder_models"),
    ]


def default_armless_cm_dir():
    """bundled default --armless-cm-dir: armless_trn{AA}_wo_{d,t,d_and_t}.cm
    files (Ozerova et al. 2024)."""
    return package_data_path("mito_cm", "armless")


# anticodon arm is the 2nd inner stem-loop (0-indexed) in a canonical cloverleaf;
# topological fact, not tunable; changing it requires a different CM.
EXPECTED_ANTICODON_STEM_INDEX = 1

# armless CM filename regex; naming follows Ozerova et al. 2024.
ARMLESS_CM_RE = re.compile(r"armless_trn(\w+)_wo_(d_and_t|d|t)\.cm$")


def index_armless_cms(cm_dir):
    """scan cm_dir for armless_trn{AA}_wo_{arm}.cm files; return {(aa, arm): path}."""
    index = _scan_cm_files(
        cm_dir, ARMLESS_CM_RE, lambda m: (m.group(1), m.group(2)), "armless"
    )
    logger.info(
        f"indexed {len(index)} armless CMs: {sorted(f'{aa}/{arm}' for aa, arm in index)}"
    )
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
    index = _scan_cm_files(
        cm_dir,
        CANONICAL_CM_RE,
        lambda m: m.group(1),
        "canonical",
        exclude=ARMLESS_CM_RE,
        warn_on_conflict=True,
    )
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
    candidates = [
        path
        for code, path in tier.items()
        if code.rstrip("0123456789") == aa_code.rstrip("0123456789")
    ]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return _pick_by_anticodon_anchor(
        header, seq, header_to_anticodon(header), candidates
    )


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


def classify_arm_loss(
    header, aligned_seq, ss_cons, expected_anticodon_index=EXPECTED_ANTICODON_STEM_INDEX
):
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
    idx, method = find_anticodon_stem_index(
        aligned_seq, elements, anticodon, expected_index=expected_anticodon_index
    )
    per_stem = [stem_complementarity(aligned_seq, ss_cons, e) for e in elements]

    result = {
        "anticodon": anticodon,
        "n_stem_loops": n,
        "anticodon_stem_index": idx,
        "anticodon_search_method": method,
        "register_offset": None,
        "per_stem_complementarity": per_stem,
        "call": "UNRESOLVED",
        "missing_arm": None,
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
        return (
            stem["n_pairs"] < MIN_STEM_PAIRS
            or stem["n_compatible"] < MIN_COMPATIBLE_PAIRS
        )

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
            if n == 3:
                # exactly D, C, T stem-loops: the 3rd is always read as the
                # T-arm. an ordinary D-C-T cloverleaf with no variable arm
                # and a real variable arm with the T-arm actually missing
                # produce the same 3-stem-loop shape; there's no structural
                # way to tell them apart (see README Limitations).
                logger.warning(
                    f"{header}: T-arm flagged absent with exactly 3 stem-loops found; "
                    "this call can't distinguish real T-arm loss from an unmodeled "
                    "variable arm - see README's 3-stem-loop limitation"
                )
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
    start = (
        elements[idx - 1]["stem3_cols"][-1] + 1
        if idx > 0
        else max(acceptor["stem5_cols"]) + 1
    )
    end = (
        elements[idx + 1]["stem5_cols"][0]
        if idx + 1 < len(elements)
        else min(acceptor["stem3_cols"])
    )
    return start, end


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
    ungapped_positions, arm_ss = _arm_full_span_subseq_and_fold(
        aligned_seq, final_seq, elem
    )

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
            ss_list[partner] = "."  # would otherwise dangle once span is cleared
        ss_list[pos] = "."

    if overridden or cleared_external:
        old_span = "".join(final_ss[p] for p in ungapped_positions)
        logger.warning(
            f"{header}: RNAfold patch overrode cmalign's own structure: "
            f"span positions {ungapped_positions[0]}-{ungapped_positions[-1]}: "
            f"old={old_span!r} -> new={arm_ss!r}; "
            f"{len(overridden)} bracket(s) inside the span replaced "
            f"({overridden}); "
            + (
                f"{len(cleared_external)} bracket(s) outside the span also cleared "
                f"to avoid dangling ({cleared_external})"
                if cleared_external
                else "none outside the span affected"
            )
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


def resolve_armless_cm(header, seq, aa_code, missing_arm, anticodon, armless_cm_index):
    """pick the correct armless CM for this sequence and arm type. if multiple
    isoacceptor CMs share the same base aa code (Leu1/Leu2, Ser1/Ser2), disambiguate
    by anticodon anchor (_pick_by_anticodon_anchor). returns None if no CM exists."""
    if aa_code is None:
        return None
    candidates = [
        path
        for (code, arm), path in armless_cm_index.items()
        if arm == missing_arm
        and code.rstrip("0123456789") == aa_code.rstrip("0123456789")
    ]
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]
    return _pick_by_anticodon_anchor(header, seq, anticodon, candidates)


def _routing_result(
    final_alignment, cm_used, diagnosis, rerouted=False, threading_failure_elem=None
):
    """assemble the dict select_cm_and_align returns at each of its exit points,
    so the shape is defined once instead of copy-pasted per branch."""
    return {
        "final_alignment": final_alignment,
        "cm_used": cm_used,
        "diagnosis": diagnosis,
        "rerouted": rerouted,
        "threading_failure_elem": threading_failure_elem,
    }


def select_cm_and_align(header, seq, canonical_cm_tiers, armless_cm_index):
    """Top-level CM selection for one mt-tRNA sequence.

    1. Align against every canonical CM tier (e.g. bacterial whole-family CM,
       then a metazoan per-AA directory). Never by raw alignment score
       (module docstring section 1).
       - Among tiers that anchor the anticodon, pick the one with the
         highest total base-pairing evidence summed across all stems
         (per_stem_complementarity's n_pairs). On a tie, prefer whichever
         tier anchored by directly matching the anticodon in one loop over
         one that only anchored by assuming its position, since a tier
         modeling an extra, here-empty stem-loop can otherwise tie a
         cleaner tier's total_pairs. A further tie (both anchored by
         position) keeps the earlier tier - not yet seen in practice.
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

    canonical_alignment = canonical_cm = diagnosis = best_key = None
    first_alignment = first_cm = first_diag = (
        None  # ultimate fallback: no tier ever anchors
    )
    for tier in canonical_cm_tiers:
        path = _resolve_canonical_for_tier(header, seq, tier)
        if path is None:
            logger.info(
                f"{header}: skipping a canonical CM tier: no CM for this amino acid there"
            )
            continue
        aln = cmalign_one(header, seq, path)
        if aln is None:
            logger.info(f"{header}: skipping a canonical CM tier: alignment failed")
            continue
        diag = classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
        if first_alignment is None:  # ultimate fallback if nothing ever anchors
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
        direct_match = diag["anticodon_search_method"] == "unique_loop_match"
        key = (total_pairs, direct_match)
        logger.debug(
            f"{header}: tier {path} anchors anticodon, total stem pairs={total_pairs}"
        )
        if canonical_alignment is None or key > best_key:
            canonical_alignment, canonical_cm, diagnosis, best_key = (
                aln,
                path,
                diag,
                key,
            )

    if canonical_alignment is None:
        canonical_alignment, canonical_cm, diagnosis = (
            first_alignment,
            first_cm,
            first_diag,
        )
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
        if elem and arm_span_has_enough_sequence(
            canonical_alignment["aligned_seq"], elem
        ):
            final_seq, _ = finalize_structure(canonical_alignment)
            if arm_is_threading_failure(
                canonical_alignment["aligned_seq"], final_seq, elem
            ):
                widened = _widen_arm_span(canonical_alignment["ss_cons"], elements, idx)
                wide_elem = dict(elem, span=widened)
                logger.info(
                    f"{header}: CM diagnosed {missing_arm}-arm missing against {canonical_cm} "
                    f"({diagnosis['call']}) but the span folds as a real hairpin "
                    f"(CM threading failure, not real arm loss); patching via RNAfold\n"
                    f"  aligned_seq={canonical_alignment['aligned_seq']}\n"
                    f"  ss_cons={canonical_alignment['ss_cons']}"
                )
                return _routing_result(
                    canonical_alignment,
                    canonical_cm,
                    diagnosis,
                    threading_failure_elem=wide_elem,
                )

    # arm loss: reroute
    logger.info(
        f"{header}: {missing_arm}-arm missing against {canonical_cm} "
        f"({diagnosis['call']}), looking for an armless CM to reroute to\n"
        f"  aligned_seq={canonical_alignment['aligned_seq']}\n"
        f"  ss_cons={canonical_alignment['ss_cons']}"
    )
    aa_code = aa_field_to_cm_code(header_to_aa(header), armless_cm_index.keys())
    anticodon = header_to_anticodon(header)
    armless_path = resolve_armless_cm(
        header, seq, aa_code, missing_arm, anticodon, armless_cm_index
    )

    if armless_path is None:
        logger.warning(
            f"{header}: {missing_arm}-arm missing ({diagnosis['call']}) "
            f"but no armless CM for aa_code={aa_code!r}; using canonical"
        )
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    armless_alignment = cmalign_one(header, seq, armless_path)
    if armless_alignment is None:
        logger.warning(
            f"{header}: armless CM realignment failed ({armless_path}); "
            f"falling back to canonical despite {missing_arm}-arm loss"
        )
        return _routing_result(canonical_alignment, canonical_cm, diagnosis)

    return _routing_result(armless_alignment, armless_path, diagnosis, rerouted=True)


# --- per-sequence worker: bundled for multiprocessing.Pool.map ---


def process_mito_record(args):
    """worker for one (header, seq) FASTA record, mito path.

    - takes a single tuple for Pool.map compatibility.
    - canonical_cm_tiers and armless_cm_index are inside that tuple: each
      worker is a fresh process, and module-level globals aren't reliably
      shared across fork vs spawn.
    - per-tier canonical CM resolution (which .cm path applies to this aa,
      if any) happens inside select_cm_and_align; see
      _resolve_canonical_for_tier."""
    header, seq, canonical_cm_tiers, armless_cm_index, debug, wc = args
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
            header,
            alignment["aligned_seq"],
            final_seq,
            final_ss,
            routing["threading_failure_elem"],
        )

    if len(final_seq) != len(seq):
        logger.warning(
            f"{header}: ungapped length {len(final_seq)} != input {len(seq)}, skipped"
        )
        return {"header": header, "rows": [], "summary": "LENGTH_MISMATCH"}

    anticodon = header_to_anticodon(header)
    if anticodon is None:
        logger.warning(f"{header}: no anticodon in header; C-stem location unreliable")

    diagnosis = routing["diagnosis"] or {}
    if wc:
        final_ss = slide_stems_to_improve_pairing(
            final_seq, final_ss, anticodon, diagnosis.get("missing_arm"), header
        )

    sprinzl = sprinzl_map(final_ss, final_seq, anticodon, diagnosis.get("missing_arm"))

    unlabeled = [i for i in range(len(final_seq)) if i not in sprinzl]
    if unlabeled:
        logger.warning(
            f"{header}: {len(unlabeled)} position(s) left without a Sprinzl "
            f"number at seq index {unlabeled}; output rows for them are blank"
        )

    cm_name = routing["cm_used"] or "NONE"
    if cm_name not in ("RNAfold", "NONE"):
        cm_name = os.path.basename(cm_name)
    summary = f"CM:{cm_name}" + (" [rerouted]" if routing["rerouted"] else "")

    logger.debug(f"{header}")
    logger.debug(f"  seq ({len(final_seq)}nt): {final_seq}")
    logger.debug(f"  ss  ({len(final_ss)}nt):  {final_ss}")
    logger.debug(
        f"  arm-loss: {diagnosis.get('call')}  "
        f"(anchor:{diagnosis.get('anticodon_search_method')}, "
        f"offset={diagnosis.get('register_offset')})"
    )
    for i, stem in enumerate(diagnosis.get("per_stem_complementarity", [])):
        logger.debug(
            f"    stem[{i}]: n_pairs={stem['n_pairs']} "
            f"n_compatible={stem['n_compatible']}"
        )
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
        rows.append(
            {
                "seq_id": header,
                "seq_index": i,
                "nucleotide": base,
                "sprinzl_position": label,
                "region": SPRINZL_REGION.get(region_key, ""),
                "cm_used": cm_name,
                "rerouted": routing["rerouted"],
                "arm_loss_call": diagnosis.get("call"),
                # dot-bracket symbol at this position; carries final_ss into the
                # TSV so scripts/visualize_ss.py can rebuild structure per record
                # without re-running cmalign (see module docstring, section 4).
                "structure": final_ss[i],
                # pre-patch CM structure and naive whole-sequence RNAfold structure,
                # same indices as final_ss (patching replaces characters in place,
                # never changes length); blank when this record wasn't a threading
                # failure, matching cm_only_ss/rnafold_only_ss being None there.
                "cm_only_structure": cm_only_ss[i] if cm_only_ss else "",
                "rnafold_only_structure": rnafold_only_ss[i] if rnafold_only_ss else "",
            }
        )

    logger.info(f"{header}: {summary}  [{diagnosis.get('call')}]")
    return {
        "header": header,
        "rows": rows,
        "summary": summary,
        "seq": final_seq,
        "ss": final_ss,
        "sprinzl": sprinzl,
        "cm_only_ss": cm_only_ss,
        "rnafold_only_ss": rnafold_only_ss,
    }


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
