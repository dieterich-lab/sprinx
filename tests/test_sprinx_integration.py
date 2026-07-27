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
  SPRINX_CANONICAL_CM=data/mito/TRNAinf-euk.cm SPRINX_ARMLESS_CM_DIR=src/sprinx/data/mito_cm/armless/ \\
  pytest test_sprinx_integration.py -v
"""
import os
import re
import shutil

import pytest

from sprinx import common, cyto, mito

MITO_CANONICAL_CM   = os.environ.get("SPRINX_CANONICAL_CM")
MITO_ARMLESS_CM_DIR = os.environ.get("SPRINX_ARMLESS_CM_DIR")
CMALIGN_OK     = shutil.which("cmalign") is not None

need_mito_cmalign = pytest.mark.skipif(
    not CMALIGN_OK or not MITO_CANONICAL_CM,
    reason="requires: cmalign in PATH, SPRINX_CANONICAL_CM env var")
need_mito_armless = pytest.mark.skipif(
    not MITO_ARMLESS_CM_DIR, reason="requires: SPRINX_ARMLESS_CM_DIR env var")
# cyto CM databases are bundled package data (no env var needed), so cyto
# tests only need cmalign itself.
need_cmalign_only = pytest.mark.skipif(not CMALIGN_OK, reason="requires cmalign in PATH")

MITO_BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "test_data_bundle.txt")
MITO_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "mito")
MITO_CM_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "sprinx", "data", "mito_cm")
MITO_BACT_CM = os.path.join(MITO_CM_DATA_DIR, "canonical", "TRNAinf-bact.cm")
MITO_METAZOA_Y_CM = os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models", "Metazoa_Y.cm")
MITO_BUNDLED_ARMLESS_CM_DIR = os.path.join(MITO_CM_DATA_DIR, "armless")
MITO_SPOMBE_FA = os.path.join(MITO_DATA_DIR, "spombe_mt.no_linker.fa")
CYTO_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "cyto")


def _load_mito_bundle_fa(key):
    text = open(MITO_BUNDLE_PATH, encoding="utf-8").read()
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
MITO_CANONICAL_SEQS = {
    "canonical_T_human.fa": "mtdbD00063518|Thr|UGU|Homo_sapiens",
    "canonical_E_human.fa": "mtdbD00063517|Glu|UUC|Homo_sapiens",
    "D_armless_human.fa":   "mtdbD00063515|Ser1|GCU|Homo_sapiens",
}


@need_mito_cmalign
def test_cmalign_one_and_finalize_real_path():
    """cmalign_one output must have equal aligned_seq/ss_cons lengths and only
    uppercase/lowercase/'-' (no '.', a multi-seq-only symbol). finalize_structure
    must then return a gap-free, balanced structure whose length matches the RNA."""
    for fa_key, seq_key in MITO_CANONICAL_SEQS.items():
        seq = _load_mito_bundle_fa(fa_key)[seq_key]
        aln = common.cmalign_one(seq_key, seq, MITO_CANONICAL_CM)
        assert aln is not None, seq_key
        assert len(aln["aligned_seq"]) == len(aln["ss_cons"]), seq_key
        assert set(aln["aligned_seq"]) <= set(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-"), seq_key
        final_seq, final_ss = common.finalize_structure(aln)
        assert len(final_seq) == len(seq) == len(final_ss), seq_key
        assert "." not in final_seq and "-" not in final_seq, seq_key
        assert final_ss.count("(") == final_ss.count(")"), seq_key


@pytest.fixture(scope="module")
def val_real_alignment():
    if not CMALIGN_OK or not MITO_CANONICAL_CM:
        pytest.skip("requires cmalign in PATH and SPRINX_CANONICAL_CM")
    seqs = _load_mito_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    return val_key, seqs[val_key], common.cmalign_one(val_key, seqs[val_key], MITO_CANONICAL_CM)


@need_mito_cmalign
def test_val_threading_failure_real_alignment(val_real_alignment):
    """the case the unit proxy could not catch: on the real single-seq span, Val's
    T-arm must be diagnosed T_OR_VAR_ARM_MISSING, pass the span check (threading
    failure, not real loss), and patch to a balanced same-length structure."""
    val_key, _seq, aln = val_real_alignment
    assert aln is not None
    elems = common.get_stem_loop_elements(aln["ss_cons"])
    t_elem = elems[-1]
    assert "T_OR_VAR_ARM_MISSING" in mito.classify_arm_loss(
        val_key, aln["aligned_seq"], aln["ss_cons"])["call"]
    assert mito.arm_span_has_enough_sequence(aln["aligned_seq"], t_elem)
    final_seq, final_ss = common.finalize_structure(aln)
    patched = mito.patch_threading_failure_arm(val_key, aln["aligned_seq"], final_seq, final_ss, t_elem)
    assert patched.count("(") == patched.count(")") and len(patched) == len(final_seq)


@need_mito_cmalign
@need_mito_armless
def test_process_one_record_populates_rnafold_only_ss_for_patched_sequences():
    """every RNAfold-patched record also carries a naive whole-sequence MFE
    fold (rnafold_only_ss, for the visualize_ss.py _RNAfoldOnly comparison
    plot), distinct from cm_only_ss (the pre-patch CM structure). neither is
    ever used for the actual patch itself, only for visual comparison."""
    seqs = _load_mito_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    armless = mito.index_armless_cms(MITO_ARMLESS_CM_DIR)
    result = mito.process_mito_record((val_key, seqs[val_key], MITO_CANONICAL_CM, armless, False))
    assert result["cm_only_ss"] is not None
    assert result["rnafold_only_ss"] is not None
    assert result["rnafold_only_ss"] != result["cm_only_ss"]
    assert len(result["rnafold_only_ss"]) == len(result["seq"])
    assert result["rnafold_only_ss"].count("(") == result["rnafold_only_ss"].count(")")


need_mito_bact_cm = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(MITO_BACT_CM)), reason="requires: cmalign, TRNAinf-bact.cm")


@need_mito_bact_cm
def test_patch_overrides_weak_pre_existing_pair_real_data():
    """a real threading-failure span (S. pombe mt-Cys's D-arm under
    TRNAinf-bact.cm, from data/mito/spombe_mt.no_linker.fa) can thread so
    weakly that only one pair survives (below MIN_STEM_PAIRS). RNAfold's fold
    of the same span agrees with that pair and extends it to a full 3bp
    D-stem: the patch applies over the pre-existing single pair rather than
    aborting, and the resulting structure leaves zero positions unlabeled."""
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(MITO_SPOMBE_FA)[header]
    aln = common.cmalign_one(header, seq, MITO_BACT_CM)
    assert aln is not None
    diag = mito.classify_arm_loss(header, aln["aligned_seq"], aln["ss_cons"])
    assert diag["missing_arm"] == "d"
    elements = common.get_stem_loop_elements(aln["ss_cons"])
    d_elem = elements[diag["anticodon_stem_index"] - 1]
    assert mito.arm_is_threading_failure(aln["aligned_seq"],
                                            common.finalize_structure(aln)[0], d_elem)
    final_seq, final_ss = common.finalize_structure(aln)
    assert final_ss.count("(") >= 1   # cmalign's own weak pair survives pre-patch
    patched = mito.patch_threading_failure_arm(header, aln["aligned_seq"], final_seq, final_ss, d_elem)
    assert patched != final_ss
    assert patched.count("(") > final_ss.count("(")
    assert patched.count("(") == patched.count(")")


@need_mito_bact_cm
def test_d_arm_patch_widens_to_recover_full_stem():
    """S. pombe mt-Cys's D-arm, once
    confirmed a threading failure, folds over the widened inter-stem domain
    (see _widen_arm_span) rather than elem['span'] alone: the narrow span
    recovers only 3bp, leaving the AD-linker 'UU' unpaired despite being
    complementary to the DC-linker's 'AA'; the wider fold recovers the full
    5bp stem, so the AD-linker ends up empty."""
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(MITO_SPOMBE_FA)[header]
    routing = mito.select_cm_and_align(header, seq, MITO_BACT_CM, {})
    assert routing["threading_failure_elem"] is not None
    aln = routing["final_alignment"]
    final_seq, final_ss = common.finalize_structure(aln)
    patched = mito.patch_threading_failure_arm(
        header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
    topo = common.parse_topology(patched)
    arms = common.locate_anticodon_stem(topo, patched, final_seq, "GCA",
                                         routing["diagnosis"]["missing_arm"])
    assert arms["linker_5"] == []
    sprinzl = common.sprinzl_map(patched, final_seq, "GCA", routing["diagnosis"]["missing_arm"])
    assert [i for i in range(len(final_seq)) if i not in sprinzl] == []


@need_mito_cmalign
@need_mito_armless
def test_select_cm_and_align_routing_and_no_unlabeled():
    """end-to-end routing on the real path, plus the no-unlabeled invariant through
    finalize + patch + sprinzl_map: canonical stays canonical, D-armless reroutes to
    a wo_d CM, Val's threading failure is patched (not rerouted)."""
    armless = mito.index_armless_cms(MITO_ARMLESS_CM_DIR)

    def _pipeline(header, seq):
        routing = mito.select_cm_and_align(header, seq, MITO_CANONICAL_CM, armless)
        aln = routing["final_alignment"]
        final_seq, final_ss = common.finalize_structure(aln)
        if routing.get("threading_failure_elem"):
            final_ss = mito.patch_threading_failure_arm(
                header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
        diag = routing["diagnosis"] or {}
        sprinzl = common.sprinzl_map(final_ss, final_seq,
                                     common.header_to_anticodon(header), diag.get("missing_arm"))
        unlabeled = [i for i in range(len(final_seq)) if i not in sprinzl]
        return routing, unlabeled

    # canonical: not rerouted
    for fa_key, seq_key in [("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
                            ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens")]:
        routing, unlabeled = _pipeline(seq_key, _load_mito_bundle_fa(fa_key)[seq_key])
        assert routing["rerouted"] is False, seq_key
        assert unlabeled == [], f"{seq_key}: unlabeled {unlabeled}"

    # D-armless: rerouted to a wo_d CM
    ser1 = "mtdbD00063515|Ser1|GCU|Homo_sapiens"
    routing, unlabeled = _pipeline(ser1, _load_mito_bundle_fa("D_armless_human.fa")[ser1])
    assert routing["rerouted"] and "wo_d" in os.path.basename(routing["cm_used"])
    assert unlabeled == []

    # Val: threading failure patched, not rerouted
    seqs = _load_mito_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    routing, unlabeled = _pipeline(val_key, seqs[val_key])
    assert routing["rerouted"] is False and routing["threading_failure_elem"] is not None
    assert unlabeled == []


@need_mito_cmalign
@need_mito_armless
def test_doubly_armless_routes_to_d_and_t_cm():
    """R. culicivorax mt-Ile (both arms absent) classifies as BOTH_ARMS_MISSING and
    routes to a d_and_t CM when one exists."""
    seqs = _load_mito_bundle_fa("both_armless_mature.fa")
    seq = seqs[next(k for k in seqs if "culicivorax" in k or "Romanomermis" in k)]
    # bundle uses id|taxon|aa|anticodon; reformat so aa_field_to_cm_code resolves.
    header = "NC_008640.1:3203-3266|Ile|GAU|Romanomermis_culicivorax"
    armless = mito.index_armless_cms(MITO_ARMLESS_CM_DIR)
    routing = mito.select_cm_and_align(header, seq, MITO_CANONICAL_CM, armless)
    diag = routing["diagnosis"]
    assert "BOTH_ARMS_MISSING" in diag["call"] or diag["missing_arm"] in ("d_and_t", "ambiguous")
    if any(arm == "d_and_t" for _, arm in armless):
        assert routing["rerouted"] and "d_and_t" in os.path.basename(routing["cm_used"])


need_mito_tiered = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(MITO_METAZOA_Y_CM) and os.path.exists(MITO_BACT_CM)),
    reason="requires: cmalign, Metazoa_Y.cm, TRNAinf-bact.cm")


@need_mito_tiered
@need_mito_armless
def test_doubly_armless_d_arm_with_zero_compatible_pairs_is_absent():
    """R. culicivorax mt-Ile's D-arm, aligned against TRNAinf-bact.cm, has 3
    non-gap column pairs but none of them are WC/wobble pairs: coincidental
    residues sitting opposite each other, not a real stem. absent() catches
    this via MIN_COMPATIBLE_PAIRS even though the raw pair count alone clears
    MIN_STEM_PAIRS, so the sequence reroutes to the d_and_t armless CM rather
    than t-only."""
    tier_dir = os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models")
    fasta = os.path.join(MITO_DATA_DIR, "both_armless.fa")
    header = "NC_008640.1:3214-3260|Ile|GAU|Romanomermis_culicivorax"
    seq = _load_fasta_file(fasta)[header]
    armless = mito.index_armless_cms(MITO_ARMLESS_CM_DIR)

    routing = mito.select_cm_and_align(
        header, seq, [MITO_BACT_CM, mito.index_canonical_cms(tier_dir)], armless)
    diag = routing["diagnosis"]
    d_arm = diag["per_stem_complementarity"][0]
    assert d_arm["n_pairs"] >= mito.MIN_STEM_PAIRS
    assert d_arm["n_compatible"] == 0
    assert diag["missing_arm"] == "d_and_t"
    assert routing["rerouted"] and "wo_d_and_t" in os.path.basename(routing["cm_used"])


need_mito_bact_and_metazoa_c = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(MITO_BACT_CM)
         and os.path.exists(os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models", "Metazoa_C.cm"))),
    reason="requires: cmalign, TRNAinf-bact.cm, Metazoa_C.cm")


@need_mito_bact_and_metazoa_c
def test_tier_prefers_fuller_anticodon_stem_thread_over_first_anchor():
    """S. pombe mt-Cys anchors cleanly against TRNAinf-bact.cm, but that CM
    only threads 3 of the anticodon stem's 5 canonical pairs; Metazoa_C.cm
    threads all 5 for the identical sequence. a short thread doesn't
    disqualify a tier outright (a real anticodon stem can be
    shorter than 5bp), but a later tier that reaches the full canonical
    count wins over one that doesn't, since accepting the short thread would
    otherwise shift the anticodon (verified via the no-unlabeled /
    anticodon-at-34-36 invariant)."""
    tier_dir = os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models")
    header = "mt-tRNA-Cys-GCA-1-1"
    seq = _load_fasta_file(MITO_SPOMBE_FA)[header]

    bact_only = mito.select_cm_and_align(header, seq, MITO_BACT_CM, {})
    idx = bact_only["diagnosis"]["anticodon_stem_index"]
    assert idx is not None
    assert bact_only["diagnosis"]["per_stem_complementarity"][idx]["n_pairs"] < mito.ANTICODON_STEM_PAIRS

    routing = mito.select_cm_and_align(header, seq, [MITO_BACT_CM, mito.index_canonical_cms(tier_dir)], {})
    assert os.path.basename(routing["cm_used"]) == "Metazoa_C.cm"
    idx2 = routing["diagnosis"]["anticodon_stem_index"]
    assert routing["diagnosis"]["per_stem_complementarity"][idx2]["n_pairs"] == mito.ANTICODON_STEM_PAIRS

    aln = routing["final_alignment"]
    final_seq, final_ss = common.finalize_structure(aln)
    if routing.get("threading_failure_elem"):
        final_ss = mito.patch_threading_failure_arm(
            header, aln["aligned_seq"], final_seq, final_ss, routing["threading_failure_elem"])
    sprinzl = common.sprinzl_map(final_ss, final_seq, "GCA", routing["diagnosis"].get("missing_arm"))
    got = "".join(final_seq[i] for i in sorted(sprinzl) if sprinzl[i] in ("34", "35", "36"))
    assert got == "GCA"
    assert [i for i in range(len(final_seq)) if i not in sprinzl] == []


@need_mito_tiered
def test_tiered_canonical_falls_back_to_bacterial():
    """S. cerevisiae mt-Tyr has a long variable loop Metazoa_Y.cm can't model
    (anchoring fails); given [Metazoa_Y, TRNAinf-bact] the bacterial tier is used."""
    header = "mtdbD00125566|Tyr|GUA|Saccharomyces cerevisiae"
    seq = ("GGAGGGAUUUUCAAUGUUGGUAGUUGGAGUUGAGCUGUAAACUCAAUGACUUAGGUCUU"
           "CAUAGGUUCAAUUCCUAUUCCCUUCA")
    # sanity: the metazoan tier alone must NOT anchor (else the fallback is moot)
    assert mito.select_cm_and_align(header, seq, MITO_METAZOA_Y_CM, {})["diagnosis"][
        "anticodon_stem_index"] is None
    routing = mito.select_cm_and_align(header, seq, [MITO_METAZOA_Y_CM, MITO_BACT_CM], {})
    assert routing["diagnosis"]["anticodon_stem_index"] is not None
    assert os.path.basename(routing["cm_used"]) == "TRNAinf-bact.cm"


need_mito_bact_armless = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(MITO_BACT_CM) and os.path.isdir(MITO_BUNDLED_ARMLESS_CM_DIR)),
    reason="requires: cmalign, TRNAinf-bact.cm, src/sprinx/data/mito_cm/armless/")


@need_mito_bact_armless
def test_all_armless_fixtures_rerouted_under_bacterial_cm():
    """TRNAinf-bact.cm models an extra variable-arm stem that broke two detection
    paths (no-shift D-arm loss; T-arm span check fooled by insert capacity). every
    ground-truth armless sequence must still reroute with the bacterial CM as the
    only canonical tier."""
    armless = mito.index_armless_cms(MITO_BUNDLED_ARMLESS_CM_DIR)
    for fa in ("D_armless.fa", "T_armless.fa"):
        seqs = _load_fasta_file(os.path.join(MITO_DATA_DIR, fa))
        not_rerouted = [h for h, s in seqs.items()
                        if not mito.select_cm_and_align(h, s, MITO_BACT_CM, armless)["rerouted"]]
        assert not_rerouted == [], f"{fa}: not rerouted: {not_rerouted}"


@need_mito_cmalign
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
    tier_dir = os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models")
    tier = mito.index_canonical_cms(tier_dir)
    assert {"L1", "L2"} <= set(tier)

    seqs = _load_mito_bundle_fa("canonical.fa")
    leu1 = next(k for k in seqs if "Leu1|UAG|Homo" in k)
    leu2 = next(k for k in seqs if "Leu2|UAA|Homo" in k)

    # simulate a GtRNAdb-style header: bare 'Leu' aa field, no isoacceptor digit.
    header1 = f"mt-tRNA-Leu-{common.header_to_anticodon(leu1)}-1-1"
    header2 = f"mt-tRNA-Leu-{common.header_to_anticodon(leu2)}-2-1"

    for header, seq in [(header1, seqs[leu1]), (header2, seqs[leu2])]:
        paths = {mito._resolve_canonical_for_tier(header, seq, tier) for _ in range(3)}
        assert len(paths) == 1, f"{header}: non-deterministic pick {paths}"
        path = paths.pop()
        assert path in tier.values()


# cyto: bundled per-isotype CM databases, no arm-loss handling.

# GtRNAdb headers with isotype-numbered aa codes (Ile2, iMet, fMet) - the
# case HEADER_TRNA_NAME_RE's 3-letter-only group used to fail to parse.
CYTO_REAL_ISOTYPE_CASES = {
    "euk": ["Homo_sapiens_tRNA-iMet-CAT-1-1", "Drosophila_melanogaster_tRNA-SeC-TCA-1-1"],
    "arch": ["Methanosarcina_barkeri_str_Fusaro_tRNA-Ile2-CAT-1-1",
             "Methanosarcina_barkeri_str_Fusaro_tRNA-iMet-CAT-1-1"],
    "bact": ["Bacillus_subtilis_subsp_subtilis_str_168_tRNA-Ile2-CAT-1-1",
             "Bacillus_subtilis_subsp_subtilis_str_168_tRNA-fMet-CAT-1-1"],
}


@need_cmalign_only
@pytest.mark.parametrize("domain", ["euk", "arch", "bact"])
def test_process_cyto_record_synthetic_consensus(domain):
    """every bundled synthetic-consensus sequence aligns against its own
    isotype CM with no unlabeled positions and the expected cm_used."""
    seqs = _load_fasta_file(os.path.join(CYTO_DATA_DIR, f"{domain}.fa"))
    cm_db_path = cyto.default_cm_db_path(domain)
    isotype_index = cyto.index_isotype_cms(cm_db_path)
    for header, seq in seqs.items():
        aa = common.header_to_aa(header)
        result = cyto.process_cyto_record((header, seq, cm_db_path, isotype_index, False))
        assert result["summary"] == f"CM:{domain}-{aa}", header
        assert len(result["rows"]) == len(seq), header
        assert all(row["sprinzl_position"] for row in result["rows"]), header


@need_cmalign_only
@pytest.mark.parametrize("domain", ["euk", "arch", "bact"])
def test_process_cyto_record_real_isotype_numbered_headers(domain):
    """real GtRNAdb headers naming Ile2/iMet/fMet/SeC align successfully -
    regression test for the HEADER_TRNA_NAME_RE fix in common.py."""
    seqs = _load_fasta_file(os.path.join(CYTO_DATA_DIR, f"{domain}_gtrnadb.fa"))
    cm_db_path = cyto.default_cm_db_path(domain)
    isotype_index = cyto.index_isotype_cms(cm_db_path)
    for header_substr in CYTO_REAL_ISOTYPE_CASES[domain]:
        header = next(h for h in seqs if header_substr in h)
        result = cyto.process_cyto_record((header, seqs[header], cm_db_path, isotype_index, False))
        assert result["summary"].startswith("CM:"), header
        assert result["rows"], header


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
