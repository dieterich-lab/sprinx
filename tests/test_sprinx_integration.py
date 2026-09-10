"""test_sprinx_integration.py - tests needing cmalign and CM files on disk

input:   tests/data/test_data_bundle.txt, data/mito/, data/cyto/, and the
         CM databases bundled under src/sprinx/data/
output:  pytest results; a test skips when its requirement is missing
usage:   SPRINX_CANONICAL_CM=data/mito/TRNAinf-euk.cm \\
             pytest tests/test_sprinx_integration.py -v
env:     cmalign in PATH, SPRINX_CANONICAL_CM, and SPRINX_ARMLESS_CM_DIR
         for the rerouting tests
notes:   single-seq cmalign differs from the multi-seq bundle proxy in
         test_sprinx_unit.py: no '.' in aligned_seq, different element
         spans, and finalize_structure length matching the RNA

The patch_threading_failure_arm crash reached production because the
multi-seq proxy could not reproduce those span differences.
"""
import multiprocessing
import os
import re
import shutil
import sys
import warnings

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

MITO_BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "data", "test_data_bundle.txt")
MITO_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "mito")
MITO_CM_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "sprinx", "data", "mito_cm")
MITO_BACT_CM = os.path.join(MITO_CM_DATA_DIR, "canonical", "TRNAinf-bact.cm")
MITO_METAZOA_Y_CM = os.path.join(MITO_CM_DATA_DIR, "canonical", "mitofinder_models", "Metazoa_Y.cm")
MITO_BUNDLED_ARMLESS_CM_DIR = os.path.join(MITO_CM_DATA_DIR, "armless")
MITO_SPOMBE_FA = os.path.join(MITO_DATA_DIR, "spombe_mt.no_linker.fa")
CYTO_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "cyto")


def _load_mito_bundle_fa(key):
    with open(MITO_BUNDLE_PATH, encoding="utf-8") as handle:
        text = handle.read()
    chunks = re.split(r"^==> (.+?) <==\n", text, flags=re.MULTILINE)[1:]
    bundle = dict(zip(chunks[0::2], chunks[1::2]))
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


# fixtures shared by the single-seq cmalign tests
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
    """on the single-seq span, Val's T-arm must be diagnosed
    T_OR_VAR_ARM_MISSING, pass the span check (threading failure rather than
    arm loss), and patch to a balanced same-length structure. The unit proxy
    cannot reach this case."""
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
    result = mito.process_mito_record(
        (val_key, seqs[val_key], MITO_CANONICAL_CM, armless, False,
         common.StructureCorrections(max_slide=0)))
    assert result["cm_only_ss"] is not None
    assert result["rnafold_only_ss"] is not None
    assert result["rnafold_only_ss"] != result["cm_only_ss"]
    assert len(result["rnafold_only_ss"]) == len(result["seq"])
    assert result["rnafold_only_ss"].count("(") == result["rnafold_only_ss"].count(")")


need_mito_bact_cm = pytest.mark.skipif(
    not (CMALIGN_OK and os.path.exists(MITO_BACT_CM)), reason="requires: cmalign, TRNAinf-bact.cm")


@need_mito_bact_cm
def test_patch_overrides_weak_pre_existing_pair_real_data():
    """a threading-failure span (S. pombe mt-Cys's D-arm under TRNAinf-bact.cm,
    from data/mito/spombe_mt.no_linker.fa) can thread so weakly that only one
    pair survives, below MIN_STEM_PAIRS. RNAfold's fold of the same span agrees
    with that pair and extends it to a full 3bp D-stem: the patch applies over
    the pre-existing single pair rather than aborting, and the result leaves
    zero positions unlabeled."""
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
    """S. pombe mt-Cys's D-arm, once confirmed a threading failure, folds over
    the widened inter-stem domain (see _widen_arm_span) rather than
    elem['span'] alone. The narrow span recovers only 3bp and leaves the
    AD-linker 'UU' unpaired despite its complementarity to the DC-linker's
    'AA'. The wider fold recovers the full 5bp stem, which empties the
    AD-linker."""
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
    """end-to-end routing on the single-seq path, plus the no-unlabeled
    invariant through finalize + patch + sprinzl_map: canonical stays canonical,
    D-armless reroutes to a wo_d CM, Val's threading failure is patched rather
    than rerouted."""
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

    # a canonical tRNA stays on its canonical CM
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
    non-gap column pairs, none of them WC/wobble: coincidental residues
    opposite each other rather than a stem. The raw pair count alone clears
    MIN_STEM_PAIRS; MIN_COMPATIBLE_PAIRS is what makes absent() catch it. The
    sequence then reroutes to the d_and_t armless CM instead of t-only."""
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
    threads only 3 of the anticodon stem's 5 canonical pairs; Metazoa_C.cm
    threads all 5 for the identical sequence. A short thread does not
    disqualify a tier outright, since an anticodon stem can be shorter than
    5bp, but a later tier reaching the full canonical count wins over one that
    does not. Accepting the short thread would shift the anticodon, checked
    here by the no-unlabeled and anticodon-at-34-36 invariants."""
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
    """for this S. cerevisiae mt-Tyr sequence, Metazoa_Y.cm aligns the region
    between the anticodon arm and T-arm as a single unpaired run, with no
    distinct stem-loop there. TRNAinf-bact.cm threads that same region as a
    paired stem with WC/wobble complementarity. Given [Metazoa_Y,
    TRNAinf-bact], the bacterial tier wins on total stem evidence. Both tiers
    anchor the anticodon by the positional tiebreak in
    find_anticodon_stem_index, which leaves stem completeness as the only
    thing separating them."""
    header = "mtdbD00125566|Tyr|GUA|Saccharomyces cerevisiae"
    seq = ("GGAGGGAUUUUCAAUGUUGGUAGUUGGAGUUGAGCUGUAAACUCAAUGACUUAGGUCUU"
           "CAUAGGUUCAAUUCCUAUUCCCUUCA")

    metazoa_only = mito.select_cm_and_align(header, seq, MITO_METAZOA_Y_CM, {})
    metazoa_elements = common.get_stem_loop_elements(metazoa_only["final_alignment"]["ss_cons"])
    assert len(metazoa_elements) == 3, "D-arm, anticodon-arm, T-arm only: no separate V-arm stem-loop"

    bact_only = mito.select_cm_and_align(header, seq, MITO_BACT_CM, {})
    bact_elements = common.get_stem_loop_elements(bact_only["final_alignment"]["ss_cons"])
    assert len(bact_elements) == 4, "D-arm, anticodon-arm, V-arm, T-arm: bact.cm models a distinct V-arm stem"
    v_arm = common.stem_complementarity(
        bact_only["final_alignment"]["aligned_seq"], bact_only["final_alignment"]["ss_cons"], bact_elements[2])
    assert v_arm["n_pairs"] > 0 and v_arm["n_compatible"] > 0

    routing = mito.select_cm_and_align(header, seq, [MITO_METAZOA_Y_CM, MITO_BACT_CM], {})
    assert routing["diagnosis"]["anticodon_stem_index"] is not None
    assert os.path.basename(routing["cm_used"]) == "TRNAinf-bact.cm"


@need_mito_tiered
def test_positional_tiebreak_yields_complete_canonical_labeling():
    """this N. tabacum plastid Thr sequence's anticodon (GGU) also occurs, by
    coincidence, in its D-loop. Against TRNAinf-bact.cm that gives 2 loop
    matches, which find_anticodon_stem_index cannot resolve on content alone.
    With 4 stem-loops present (D, anticodon, V, T), index 1 is unambiguous
    whatever the loop content, and the anticodon anchors there. The downstream
    Sprinzl assignment then comes out complete: a 7bp acceptor stem, every
    position labeled, anticodon at 34-36."""
    header = "NtK326-Chl-tRNA-Thr-GGT-1-1-NC_001879:33186-33257 (+) 72 bp Sc:58.2 -"
    seq = "GCCCUUUUAACUCAGUGGUAGAGUAACGCCAUGGUAAGGCGUAAGUCAUCGGUUCAAAUCCGAUAAGGGGCU"

    aln = common.cmalign_one(header, seq, MITO_BACT_CM)
    elements = common.get_stem_loop_elements(aln["ss_cons"])
    assert common.find_anticodon_stem_index(aln["aligned_seq"], elements, "GGU") == (
        None, "ambiguous_2_loop_matches")

    routing = mito.select_cm_and_align(header, seq, MITO_BACT_CM, {})
    diag = routing["diagnosis"]
    assert diag["anticodon_stem_index"] == 1
    assert not routing["rerouted"]

    final_seq, final_ss = common.finalize_structure(routing["final_alignment"])
    sprinzl = common.sprinzl_map(final_ss, final_seq, "GGU", diag["missing_arm"])
    assert [i for i in range(len(final_seq)) if i not in sprinzl] == []
    assert [sprinzl[i] for i in range(7)] == ["1", "2", "3", "4", "5", "6", "7"]
    got_anticodon = "".join(final_seq[i] for i in sorted(sprinzl)
                            if sprinzl[i] in ("34", "35", "36"))
    assert got_anticodon == "GGU"


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
    """a GtRNAdb-style header never carries an isoacceptor digit: bare 'Leu'
    covers both anticodons, and a per-AA tier with separate L1/L2 CMs has to
    resolve that code some way other than a dict lookup. Metazoa_L1.cm and
    L2.cm both anchor either Leu anticodon equally well - the isoacceptor
    filenames label the two CMs without encoding a structural difference the
    anchor check could see. So the assertions here are only that resolution
    returns some candidate rather than None, and that it is deterministic: a
    flaky pick would make Sprinzl output non-reproducible."""
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

# GtRNAdb headers with isotype-numbered aa codes (Ile2, iMet, fMet), which
# HEADER_TRNA_NAME_RE's aa group must accept alongside plain 3-letter codes.
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
        result = cyto.process_cyto_record(
            (header, seq, cm_db_path, isotype_index, False,
             common.StructureCorrections(max_slide=0)))
        assert result["summary"] == f"CM:{domain}-{aa}", header
        assert len(result["rows"]) == len(seq), header
        assert all(row["sprinzl_position"] for row in result["rows"]), header


@need_cmalign_only
@pytest.mark.parametrize("domain", ["euk", "arch", "bact"])
def test_process_cyto_record_real_isotype_numbered_headers(domain):
    """GtRNAdb headers naming Ile2/iMet/fMet/SeC align successfully. Guards
    HEADER_TRNA_NAME_RE against a regression on isotype-numbered aa codes."""
    seqs = _load_fasta_file(os.path.join(CYTO_DATA_DIR, f"{domain}_gtrnadb.fa"))
    cm_db_path = cyto.default_cm_db_path(domain)
    isotype_index = cyto.index_isotype_cms(cm_db_path)
    for header_substr in CYTO_REAL_ISOTYPE_CASES[domain]:
        header = next(h for h in seqs if header_substr in h)
        result = cyto.process_cyto_record(
            (header, seqs[header], cm_db_path, isotype_index, False,
             common.StructureCorrections(max_slide=0)))
        assert result["summary"].startswith("CM:"), header
        assert result["rows"], header


CONSERVED_POSITIONS_PATH = os.path.join(os.path.dirname(__file__), "data", "conserved_positions.tsv")


def _load_conserved_positions():
    """(position, base, min_fraction) rows. min_fraction is None for a row
    marked "log". Such a row is reported, never asserted."""
    rows = []
    with open(CONSERVED_POSITIONS_PATH, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            position, base, min_fraction, _measured = line.split()
            rows.append((position, base,
                         None if min_fraction == "log" else float(min_fraction)))
    return rows


@pytest.fixture(scope="module")
def euk_gtrnadb_labels():
    """name -> (base, sprinzl_position) per sequence position in 5'->3' order,
    for the whole euk GtRNAdb set, with '' where sprinx assigned no label.
    Labeling all 531 sequences takes about 12 seconds, shared across the
    conserved-base and CM-agreement tests. The full set is deliberate: the
    conserved-base fractions are measured over natural sequence diversity, and
    collapsing the input would change what they mean."""
    seqs = _load_fasta_file(os.path.join(CYTO_DATA_DIR, "euk_gtrnadb.fa"))
    cm_db_path = cyto.default_cm_db_path("euk")
    isotype_index = cyto.index_isotype_cms(cm_db_path)
    tasks = [(header, seq, cm_db_path, isotype_index, False, common.StructureCorrections())
             for header, seq in seqs.items()]
    with multiprocessing.Pool(4) as pool:
        results = pool.map(cyto.process_cyto_record, tasks)

    return {header.split()[0]: [(row["nucleotide"], row["sprinzl_position"])
                                for row in result["rows"]]
            for header, result in zip(seqs, results)}


@pytest.fixture(scope="module")
def euk_gtrnadb_bases_by_position(euk_gtrnadb_labels):
    """sprinzl_position -> bases assigned to it across the whole euk GtRNAdb set."""
    by_position = {}
    for rows in euk_gtrnadb_labels.values():
        for base, label in rows:
            if label:
                by_position.setdefault(label, []).append(base)
    return by_position


@need_cmalign_only
def test_conserved_positions_carry_expected_bases(euk_gtrnadb_bases_by_position):
    """labels that slip off their column stop landing on the base their Sprinzl
    position is known to carry, which shows up here as the fraction dropping
    below min_fraction. Every position is reported at once, since one shift
    usually drags its neighbours with it.

    Cytosolic input only. A row marked "log" is warned, not asserted. None of
    these fractions is measured on mt-tRNA, where several of the bases are not
    conserved. None anchors an assignment."""
    failures = []
    for position, base, min_fraction in _load_conserved_positions():
        observed = euk_gtrnadb_bases_by_position.get(position, [])
        if not observed:
            failures.append(f"{position}: no sequence was labeled with this position")
            continue
        fraction = observed.count(base) / len(observed)
        if min_fraction is None:
            warnings.warn(f"conserved position {position}: {base} in {fraction:.3f} "
                          f"of {len(observed)} cytosolic sequences (reported, not gated)",
                          stacklevel=2)
            continue
        if fraction < min_fraction:
            failures.append(f"{position}: {base} in {fraction:.3f} of {len(observed)} "
                            f"sequences, under the {min_fraction} floor")
    assert not failures, "conserved positions below their floor:\n  " + "\n  ".join(failures)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
