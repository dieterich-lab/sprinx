"""
test_sprinx_integration.py: integration tests requiring cmalign and real CM files.

these exercise the real single-seq cmalign path, which differs from the multi-seq
bundle proxy used in test_sprinx_unit.py: single-seq aligned_seq has no '.' chars,
element spans differ, and finalize_structure output length matches the real RNA.
the production crash in patch_threading_failure_arm was invisible to the unit proxy
precisely because of these span differences.

requirements (all must be set to run any test here):
  cmalign in PATH
  SPRINX_CANONICAL_CM   : path to the canonical CM (e.g. TRNAinf-euk.cm)
  SPRINX_ARMLESS_CM_DIR : directory of armless CMs (only for rerouting tests)

run:
  SPRINX_CANONICAL_CM=data/TRNAinf-euk.cm SPRINX_ARMLESS_CM_DIR=data/truncated_cm/ \\
  pytest test_sprinx_integration.py -v
"""
import os
import re
import shutil

import pytest

from sprinx import label as sprinx

CANONICAL_CM   = os.environ.get("SPRINX_CANONICAL_CM")
ARMLESS_CM_DIR = os.environ.get("SPRINX_ARMLESS_CM_DIR")
CMALIGN_OK     = shutil.which("cmalign") is not None

need_cmalign = pytest.mark.skipif(
    not CMALIGN_OK or not CANONICAL_CM,
    reason="requires: cmalign in PATH, SPRINX_CANONICAL_CM env var")
need_armless = pytest.mark.skipif(
    not ARMLESS_CM_DIR, reason="requires: SPRINX_ARMLESS_CM_DIR env var")

BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "test_data_bundle.txt")
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
BACT_CM = os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE", "TRNAinf-bact.cm")
METAZOA_Y_CM = os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE", "Metazoa_Y.cm")
TRUNCATED_CM_DIR = os.path.join(DATA_DIR, "truncated_cm")
SPOMBE_MT_FA = os.path.join(DATA_DIR, "spombe_mt.no_linker.fa")


def _load_bundle_fa(key):
    text = open(BUNDLE_PATH, encoding="utf-8").read()
    chunks = re.split(r"^==> (.+?) <==\n", text, flags=re.MULTILINE)[1:]
    bundle = {name: content for name, content in zip(chunks[0::2], chunks[1::2])}
    seqs, cur = {}, None
    for line in bundle[key].splitlines():
        if line.startswith(">"):
            cur = line[1:].strip()
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip().upper().replace("T", "U")
    return seqs


def _load_fasta_file(path):
    seqs, cur = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                cur = line[1:].strip()
                seqs[cur] = ""
            elif cur:
                seqs[cur] += line.strip().upper().replace("T", "U")
    return seqs


# fixtures used across the real-path tests
CANONICAL_SEQS = {
    "canonical_T_human.fa": "mtdbD00063518|Thr|UGU|Homo_sapiens",
    "canonical_E_human.fa": "mtdbD00063517|Glu|UUC|Homo_sapiens",
    "D_armless_human.fa":   "mtdbD00063515|Ser1|GCU|Homo_sapiens",
}


@need_cmalign
def test_cmalign_one_and_finalize_real_path():
    """cmalign_one output must have equal aligned_seq/ss_cons lengths and only
    uppercase/lowercase/'-' (no '.', a multi-seq-only symbol). finalize_structure
    must then return a gap-free, balanced structure whose length matches the RNA."""
    for fa_key, seq_key in CANONICAL_SEQS.items():
        seq = _load_bundle_fa(fa_key)[seq_key]
        aln = sprinx.cmalign_one(seq_key, seq, CANONICAL_CM)
        assert aln is not None, seq_key
        assert len(aln["aligned_seq"]) == len(aln["ss_cons"]), seq_key
        assert set(aln["aligned_seq"]) <= set(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-"), seq_key
        final_seq, final_ss = sprinx.finalize_structure(aln)
        assert len(final_seq) == len(seq) == len(final_ss), seq_key
        assert "." not in final_seq and "-" not in final_seq, seq_key
        assert final_ss.count("(") == final_ss.count(")"), seq_key


@pytest.fixture(scope="module")
def val_real_alignment():
    if not CMALIGN_OK or not CANONICAL_CM:
        pytest.skip("requires cmalign in PATH and SPRINX_CANONICAL_CM")
    seqs = _load_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    return val_key, seqs[val_key], sprinx.cmalign_one(val_key, seqs[val_key], CANONICAL_CM)


@need_cmalign
def test_val_threading_failure_real_alignment(val_real_alignment):
    """the case the unit proxy could not catch: on the real single-seq span, Val's
    T-arm must be diagnosed T_OR_VAR_ARM_MISSING, pass the span check (threading
    failure, not real loss), and patch to a balanced same-length structure."""
    val_key, _seq, aln = val_real_alignment
    assert aln is not None
    elems = sprinx.get_stem_loop_elements(aln["ss_cons"])
    t_elem = elems[-1]
    assert "T_OR_VAR_ARM_MISSING" in sprinx.classify_arm_loss(
        val_key, aln["aligned_seq"], aln["ss_cons"])["call"]
    assert sprinx.arm_span_has_enough_sequence(aln["aligned_seq"], t_elem)
    final_seq, final_ss = sprinx.finalize_structure(aln)
    patched = sprinx.patch_threading_failure_arm(val_key, aln["aligned_seq"], final_seq, final_ss, t_elem)
    assert patched.count("(") == patched.count(")") and len(patched) == len(final_seq)


@need_cmalign
@need_armless
def test_process_one_record_populates_rnafold_only_ss_for_patched_sequences():
    """every RNAfold-patched record also carries a naive whole-sequence MFE
    fold (rnafold_only_ss, for the visualize_ss.py _RNAfoldOnly comparison
    plot), distinct from cm_only_ss (the pre-patch CM structure). neither is
    ever used for the actual patch itself, only for visual comparison."""
    seqs = _load_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    armless = sprinx.index_armless_cms(ARMLESS_CM_DIR)
    result = sprinx.process_one_record((val_key, seqs[val_key], CANONICAL_CM, armless, False))
    assert result["cm_only_ss"] is not None
    assert result["rnafold_only_ss"] is not None
    assert result["rnafold_only_ss"] != result["cm_only_ss"]
    assert len(result["rnafold_only_ss"]) == len(result["seq"])
    assert result["rnafold_only_ss"].count("(") == result["rnafold_only_ss"].count(")")


need_bact = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(BACT_CM)), reason="requires: cmalign, TRNAinf-bact.cm")


@need_bact
def test_patch_overrides_weak_pre_existing_pair_real_data():
    """a real threading-failure span (S. pombe mt-Cys's D-arm under
    TRNAinf-bact.cm, from data/spombe_mt.no_linker.fa) can thread so
    weakly that only one pair survives (below MIN_STEM_PAIRS). RNAfold's fold
    of the same span agrees with that pair and extends it to a full 3bp
    D-stem: the patch applies over the pre-existing single pair rather than
    aborting, and the resulting structure leaves zero positions unlabeled."""
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(SPOMBE_MT_FA)[header]
    aln = sprinx.cmalign_one(header, seq, BACT_CM)
    assert aln is not None
    diag = sprinx.classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
    assert diag["missing_arm"] == "d"
    elements = sprinx.get_stem_loop_elements(aln["ss_cons"])
    d_elem = elements[diag["anticodon_stem_index"] - 1]
    assert sprinx.arm_is_threading_failure(aln["aligned_seq"],
                                            sprinx.finalize_structure(aln)[0], d_elem)
    final_seq, final_ss = sprinx.finalize_structure(aln)
    assert final_ss.count("(") >= 1   # cmalign's own weak pair survives pre-patch
    patched = sprinx.patch_threading_failure_arm(header, aln["aligned_seq"], final_seq, final_ss, d_elem)
    assert patched != final_ss
    assert patched.count("(") > final_ss.count("(")
    assert patched.count("(") == patched.count(")")


@need_bact
def test_d_arm_patch_widens_to_recover_full_stem():
    """S. pombe mt-Cys's D-arm, once
    confirmed a threading failure, folds over the widened inter-stem domain
    (see _widen_arm_span) rather than elem['span'] alone: the narrow span
    recovers only 3bp, leaving the AD-linker 'UU' unpaired despite being
    complementary to the DC-linker's 'AA'; the wider fold recovers the full
    5bp stem, so the AD-linker ends up empty."""
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(SPOMBE_MT_FA)[header]
    routing = sprinx.select_cm_and_align(header, seq, BACT_CM, {})
    assert routing["threading_failure_elem"] is not None
    aln = routing["final_alignment"]
    final_seq, final_ss = sprinx.finalize_structure(aln)
    patched = sprinx.patch_threading_failure_arm(
        header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
    topo = sprinx.parse_topology(patched)
    arms = sprinx.locate_anticodon_stem(topo, patched, final_seq, "GCA",
                                         routing["diagnosis"]["missing_arm"])
    assert arms["linker_5"] == []
    sprinzl = sprinx.sprinzl_map(patched, final_seq, "GCA", routing["diagnosis"]["missing_arm"])
    assert [i for i in range(len(final_seq)) if i not in sprinzl] == []


@need_cmalign
@need_armless
def test_select_cm_and_align_routing_and_no_unlabeled():
    """end-to-end routing on the real path, plus the no-unlabeled invariant through
    finalize + patch + sprinzl_map: canonical stays canonical, D-armless reroutes to
    a wo_d CM, Val's threading failure is patched (not rerouted)."""
    armless = sprinx.index_armless_cms(ARMLESS_CM_DIR)

    def _pipeline(header, seq):
        routing = sprinx.select_cm_and_align(header, seq, CANONICAL_CM, armless)
        aln = routing["final_alignment"]
        final_seq, final_ss = sprinx.finalize_structure(aln)
        if routing.get("threading_failure_elem"):
            final_ss = sprinx.patch_threading_failure_arm(
                header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
        diag = routing["diagnosis"] or {}
        sprinzl = sprinx.sprinzl_map(final_ss, final_seq,
                                     sprinx.header_to_anticodon(header), diag.get("missing_arm"))
        unlabeled = [i for i in range(len(final_seq)) if i not in sprinzl]
        return routing, unlabeled

    # canonical: not rerouted
    for fa_key, seq_key in [("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
                            ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens")]:
        routing, unlabeled = _pipeline(seq_key, _load_bundle_fa(fa_key)[seq_key])
        assert routing["rerouted"] is False, seq_key
        assert unlabeled == [], f"{seq_key}: unlabeled {unlabeled}"

    # D-armless: rerouted to a wo_d CM
    ser1 = "mtdbD00063515|Ser1|GCU|Homo_sapiens"
    routing, unlabeled = _pipeline(ser1, _load_bundle_fa("D_armless_human.fa")[ser1])
    assert routing["rerouted"] and "wo_d" in os.path.basename(routing["cm_used"])
    assert unlabeled == []

    # Val: threading failure patched, not rerouted
    seqs = _load_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    routing, unlabeled = _pipeline(val_key, seqs[val_key])
    assert routing["rerouted"] is False and routing["threading_failure_elem"] is not None
    assert unlabeled == []


@need_cmalign
@need_armless
def test_doubly_armless_routes_to_d_and_t_cm():
    """R. culicivorax mt-Ile (both arms absent) classifies as BOTH_ARMS_MISSING and
    routes to a d_and_t CM when one exists."""
    seqs = _load_bundle_fa("both_armless_mature.fa")
    seq = seqs[next(k for k in seqs if "culicivorax" in k or "Romanomermis" in k)]
    # bundle uses id|taxon|aa|anticodon; reformat so aa_field_to_cm_code resolves.
    header = "NC_008640.1:3203-3266|Ile|GAU|Romanomermis_culicivorax"
    armless = sprinx.index_armless_cms(ARMLESS_CM_DIR)
    routing = sprinx.select_cm_and_align(header, seq, CANONICAL_CM, armless)
    diag = routing["diagnosis"]
    assert "BOTH_ARMS_MISSING" in diag["call"] or diag["missing_arm"] in ("d_and_t", "ambiguous")
    if any(arm == "d_and_t" for _, arm in armless):
        assert routing["rerouted"] and "d_and_t" in os.path.basename(routing["cm_used"])


need_tiered = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(METAZOA_Y_CM) and os.path.exists(BACT_CM)),
    reason="requires: cmalign, Metazoa_Y.cm, TRNAinf-bact.cm")


@need_tiered
@need_armless
def test_doubly_armless_d_arm_with_zero_compatible_pairs_is_absent():
    """R. culicivorax mt-Ile's D-arm, aligned against TRNAinf-bact.cm, has 3
    non-gap column pairs but none of them are WC/wobble pairs: coincidental
    residues sitting opposite each other, not a real stem. absent() catches
    this via MIN_COMPATIBLE_PAIRS even though the raw pair count alone clears
    MIN_STEM_PAIRS, so the sequence reroutes to the d_and_t armless CM rather
    than t-only."""
    tier_dir = os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE")
    fasta = os.path.join(DATA_DIR, "both_armless.fa")
    header = "NC_008640.1:3214-3260|Ile|GAU|Romanomermis_culicivorax"
    seq = _load_fasta_file(fasta)[header]
    armless = sprinx.index_armless_cms(ARMLESS_CM_DIR)

    routing = sprinx.select_cm_and_align(
        header, seq, [BACT_CM, sprinx.index_canonical_cms(tier_dir)], armless)
    diag = routing["diagnosis"]
    d_arm = diag["per_stem_complementarity"][0]
    assert d_arm["n_pairs"] >= sprinx.MIN_STEM_PAIRS
    assert d_arm["n_compatible"] == 0
    assert diag["missing_arm"] == "d_and_t"
    assert routing["rerouted"] and "wo_d_and_t" in os.path.basename(routing["cm_used"])


need_bact_and_metazoa_c = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(BACT_CM)
         and os.path.exists(os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE", "Metazoa_C.cm"))),
    reason="requires: cmalign, TRNAinf-bact.cm, Metazoa_C.cm")


@need_bact_and_metazoa_c
def test_tier_prefers_fuller_anticodon_stem_thread_over_first_anchor():
    """S. pombe mt-Cys anchors cleanly against TRNAinf-bact.cm, but that CM
    only threads 3 of the anticodon stem's 5 canonical pairs; Metazoa_C.cm
    threads all 5 for the identical sequence. a short thread doesn't
    disqualify a tier outright (a real anticodon stem can be
    shorter than 5bp), but a later tier that reaches the full canonical
    count wins over one that doesn't, since accepting the short thread would
    otherwise shift the anticodon (verified via the no-unlabeled /
    anticodon-at-34-36 invariant)."""
    tier_dir = os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE")
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(SPOMBE_MT_FA)[header]

    bact_only = sprinx.select_cm_and_align(header, seq, BACT_CM, {})
    idx = bact_only["diagnosis"]["anticodon_stem_index"]
    assert idx is not None
    assert bact_only["diagnosis"]["per_stem_complementarity"][idx]["n_pairs"] < sprinx.ANTICODON_STEM_PAIRS

    routing = sprinx.select_cm_and_align(header, seq, [BACT_CM, sprinx.index_canonical_cms(tier_dir)], {})
    assert os.path.basename(routing["cm_used"]) == "Metazoa_C.cm"
    idx2 = routing["diagnosis"]["anticodon_stem_index"]
    assert routing["diagnosis"]["per_stem_complementarity"][idx2]["n_pairs"] == sprinx.ANTICODON_STEM_PAIRS

    aln = routing["final_alignment"]
    final_seq, final_ss = sprinx.finalize_structure(aln)
    if routing.get("threading_failure_elem"):
        final_ss = sprinx.patch_threading_failure_arm(
            header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
    sprinzl = sprinx.sprinzl_map(final_ss, final_seq, "GCA", routing["diagnosis"].get("missing_arm"))
    got = "".join(final_seq[i] for i in sorted(sprinzl) if sprinzl[i] in ("34", "35", "36"))
    assert got == "GCA"
    assert [i for i in range(len(final_seq)) if i not in sprinzl] == []


@need_tiered
def test_tiered_canonical_falls_back_to_bacterial():
    """S. cerevisiae mt-Tyr has a long variable loop Metazoa_Y.cm can't model
    (anchoring fails); given [Metazoa_Y, TRNAinf-bact] the bacterial tier is used."""
    header = "mtdbD00125566|Tyr|GUA|Saccharomyces cerevisiae"
    seq = ("GGAGGGAUUUUCAAUGUUGGUAGUUGGAGUUGAGCUGUAAACUCAAUGACUUAGGUCUU"
           "CAUAGGUUCAAUUCCUAUUCCCUUCA")
    # sanity: the metazoan tier alone must NOT anchor (else the fallback is moot)
    assert sprinx.select_cm_and_align(header, seq, METAZOA_Y_CM, {})["diagnosis"][
        "anticodon_stem_index"] is None
    routing = sprinx.select_cm_and_align(header, seq, [METAZOA_Y_CM, BACT_CM], {})
    assert routing["diagnosis"]["anticodon_stem_index"] is not None
    assert os.path.basename(routing["cm_used"]) == "TRNAinf-bact.cm"


need_bact_armless = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(BACT_CM) and os.path.isdir(TRUNCATED_CM_DIR)),
    reason="requires: cmalign, TRNAinf-bact.cm, data/truncated_cm/")


@need_bact_armless
def test_all_armless_fixtures_rerouted_under_bacterial_cm():
    """TRNAinf-bact.cm models an extra variable-arm stem that broke two detection
    paths (no-shift D-arm loss; T-arm span check fooled by insert capacity). every
    ground-truth armless sequence must still reroute with the bacterial CM as the
    only canonical tier."""
    armless = sprinx.index_armless_cms(TRUNCATED_CM_DIR)
    for fa in ("D_armless.fa", "T_armless.fa"):
        seqs = _load_fasta_file(os.path.join(DATA_DIR, fa))
        not_rerouted = [h for h, s in seqs.items()
                        if not sprinx.select_cm_and_align(h, s, BACT_CM, armless)["rerouted"]]
        assert not_rerouted == [], f"{fa}: not rerouted: {not_rerouted}"


def test_resolve_canonical_for_tier():
    # a plain path applies to every aa; a dict resolves by aa or returns None.
    assert sprinx._resolve_canonical_for_tier("any|header|x|y", "ACGU", "/p.cm") == "/p.cm"
    tier = {"A": "/models/Ala.cm", "V": "/models/Val.cm"}
    assert sprinx._resolve_canonical_for_tier("id|Ala|UGC|taxon", "ACGU", tier) == "/models/Ala.cm"
    assert sprinx._resolve_canonical_for_tier("id|Trp|UCA|taxon", "ACGU", tier) is None


@need_cmalign
def test_resolve_canonical_for_tier_disambiguates_bare_isoacceptor_by_anticodon():
    """a GtRNAdb-style header never carries an isoacceptor digit (bare 'Leu'
    covers both anticodons), so a per-AA tier with separate L1/L2 CMs has to
    resolve a bare aa code some way other than a direct dict lookup. real
    Metazoa_L1.cm/L2.cm both structurally anchor either real Leu anticodon
    equally well, since the filename split isn't itself something the anchor
    check can see, consistent with isoacceptor filenames being arbitrary
    labels rather than a structural distinction. so this can't assert which
    specific file comes back; what matters instead is that resolution never
    fails silently (a bare code always returns *some* real candidate, not
    None) and is deterministic (same header+seq -> same CM every call, since
    a flaky pick would make Sprinzl output non-reproducible)."""
    tier_dir = os.path.join(DATA_DIR, "full_tRNAs_mitofinder_tRNAScanSE")
    tier = sprinx.index_canonical_cms(tier_dir)
    assert {"L1", "L2"} <= set(tier)

    seqs = _load_bundle_fa("canonical.fa")
    leu1 = next(k for k in seqs if "Leu1|UAG|Homo" in k)
    leu2 = next(k for k in seqs if "Leu2|UAA|Homo" in k)

    # simulate a GtRNAdb-style header: bare 'Leu' aa field, no isoacceptor digit.
    header1 = f"mt-tRNA-Leu-{sprinx.header_to_anticodon(leu1)}-1-1"
    header2 = f"mt-tRNA-Leu-{sprinx.header_to_anticodon(leu2)}-2-1"

    for header, seq in [(header1, seqs[leu1]), (header2, seqs[leu2])]:
        paths = {sprinx._resolve_canonical_for_tier(header, seq, tier) for _ in range(3)}
        assert len(paths) == 1, f"{header}: non-deterministic pick {paths}"
        path = paths.pop()
        assert path in tier.values()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
