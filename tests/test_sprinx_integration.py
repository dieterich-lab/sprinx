"""
test_sprinx_integration.py -- integration tests requiring cmalign and real CM files.

these tests exercise the actual single-seq cmalign output path, which differs from
the multi-seq bundle proxy used in test_sprinx_unit.py:
  - single-seq aligned_seq has no '.' characters (only uppercase/lowercase/'-')
  - element spans in the real gapped SS_cons differ from multi-seq spans
  - finalize_structure output length matches the actual RNA, not a hardcoded proxy

the production crash in patch_threading_failure_arm was not caught by unit tests
precisely because the multi-seq proxy produced different element spans than the real
single-seq alignment; an integration test running cmalign against the real CM would
have caught it immediately.

requirements (all must be set to run any test here):
  cmalign in PATH
  SPRINX_CANONICAL_CM   : path to the canonical CM (e.g. TRNAinf-euk.cm)
  SPRINX_ARMLESS_CM_DIR : path to directory of armless CMs (armless_trn*_wo_*.cm)
                          only required for rerouting tests; others run without it.

run:
  SPRINX_CANONICAL_CM=data/TRNAinf-euk.cm \\
  SPRINX_ARMLESS_CM_DIR=data/truncated_cm/ \\
  pytest test_sprinx_integration.py -v
"""
import os
import re
import shutil
import sys

import pytest
import RNA

sys.path.insert(0, os.path.dirname(__file__))
import sprinx  # pylint: disable=wrong-import-position

CANONICAL_CM  = os.environ.get("SPRINX_CANONICAL_CM")
ARMLESS_CM_DIR = os.environ.get("SPRINX_ARMLESS_CM_DIR")
CMALIGN_OK    = shutil.which("cmalign") is not None

# per-test skip markers
need_cmalign = pytest.mark.skipif(
    not CMALIGN_OK or not CANONICAL_CM,
    reason="requires: cmalign in PATH, SPRINX_CANONICAL_CM env var"
)
need_armless = pytest.mark.skipif(
    not ARMLESS_CM_DIR,
    reason="requires: SPRINX_ARMLESS_CM_DIR env var"
)

# -----------------------------------------------------------------------
# sequences used across tests; sourced from bundle FASTA files so they
# match what the unit tests use, without hardcoding by hand.
# -----------------------------------------------------------------------

BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "test_data_bundle.txt")


def _load_bundle_fa(bundle_key):
    """parse one FASTA block from test_data_bundle.txt into {header: seq}."""
    text = open(BUNDLE_PATH, encoding="utf-8").read()
    chunks = re.split(r"^==> (.+?) <==\n", text, flags=re.MULTILINE)[1:]
    bundle = {name: content for name, content in zip(chunks[0::2], chunks[1::2])}
    seqs = {}
    cur = None
    for line in bundle[bundle_key].splitlines():
        if line.startswith(">"):
            cur = line[1:].strip()
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip().upper().replace("T", "U")
    return seqs


# -----------------------------------------------------------------------
# cmalign_one: alignment parsing correctness
# -----------------------------------------------------------------------

class TestCmalignOneParsing:
    """verify that cmalign_one correctly parses the cmalign Stockholm output:
    aligned_seq and ss_cons lengths must match (a length mismatch here means
    finalize_structure will silently truncate one of them via zip, producing
    a wrong-length final_seq or final_ss)."""

    @need_cmalign
    @pytest.mark.parametrize("fa_key, seq_key", [
        ("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
        ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens"),
        ("canonical_L1_human.fa", "mtdbD00063516|Leu1|UAG|Homo_sapiens"),
    ])
    def test_aligned_seq_and_ss_cons_same_length(self, fa_key, seq_key):
        seqs = _load_bundle_fa(fa_key)
        seq = seqs[seq_key]
        aln = sprinx.cmalign_one(seq_key, seq, CANONICAL_CM)
        assert aln is not None, f"cmalign_one returned None for {seq_key}"
        assert len(aln["aligned_seq"]) == len(aln["ss_cons"]), (
            f"{seq_key}: aligned_seq ({len(aln['aligned_seq'])}) and "
            f"ss_cons ({len(aln['ss_cons'])}) lengths differ"
        )

    @need_cmalign
    def test_aligned_seq_character_set_single_seq(self):
        """single-seq cmalign output must use only uppercase (match), lowercase
        (insert), and '-' (deletion). no '.' characters. this is the critical
        difference from multi-seq alignments that the unit test proxy cannot check."""
        seqs = _load_bundle_fa("canonical_E_human.fa")
        seq_key = "mtdbD00063517|Glu|UUC|Homo_sapiens"
        aln = sprinx.cmalign_one(seq_key, seqs[seq_key], CANONICAL_CM)
        assert aln is not None
        valid = set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz-")
        invalid = set(aln["aligned_seq"]) - valid
        assert not invalid, (
            f"single-seq aligned_seq contains unexpected characters: {invalid}. "
            f"'.' is a multi-seq-only gap symbol and must not appear here."
        )


# -----------------------------------------------------------------------
# finalize_structure: real single-seq cmalign output
# -----------------------------------------------------------------------

class TestFinalizeStructureIntegration:
    """test finalize_structure on actual single-seq cmalign output. the unit
    test uses a multi-seq bundle proxy, which has different gap character
    distributions and element spans. these tests use the real code path."""

    @need_cmalign
    @pytest.mark.parametrize("fa_key, seq_key", [
        ("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
        ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens"),
        ("D_armless_human.fa",   "mtdbD00063515|Ser1|GCU|Homo_sapiens"),
    ])
    def test_seq_and_ss_same_length_as_input(self, fa_key, seq_key):
        seqs = _load_bundle_fa(fa_key)
        seq = seqs[seq_key]
        aln = sprinx.cmalign_one(seq_key, seq, CANONICAL_CM)
        assert aln is not None
        final_seq, final_ss = sprinx.finalize_structure(aln)
        assert len(final_seq) == len(seq), (
            f"{seq_key}: finalize_structure output length ({len(final_seq)}) "
            f"!= input length ({len(seq)})"
        )
        assert len(final_ss) == len(final_seq), (
            f"{seq_key}: final_ss ({len(final_ss)}) and final_seq ({len(final_seq)}) differ"
        )

    @need_cmalign
    @pytest.mark.parametrize("fa_key, seq_key", [
        ("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
        ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens"),
        ("D_armless_human.fa",   "mtdbD00063515|Ser1|GCU|Homo_sapiens"),
    ])
    def test_output_is_clean_and_balanced(self, fa_key, seq_key):
        """final_seq must contain no gap characters; final_ss must be balanced."""
        seqs = _load_bundle_fa(fa_key)
        seq = seqs[seq_key]
        aln = sprinx.cmalign_one(seq_key, seq, CANONICAL_CM)
        assert aln is not None
        final_seq, final_ss = sprinx.finalize_structure(aln)
        assert "." not in final_seq and "-" not in final_seq, (
            f"{seq_key}: gap chars in final_seq"
        )
        assert final_ss.count("(") == final_ss.count(")"), (
            f"{seq_key}: unbalanced final_ss ({final_ss.count('(')} opens, "
            f"{final_ss.count(')')} closes)"
        )


# -----------------------------------------------------------------------
# Val threading failure: real single-seq path
#
# this is the test that would have caught the production crash. the unit test
# uses a multi-seq proxy with different element spans; here we use the real
# single-seq cmalign output so the T-arm element span reflects production.
# -----------------------------------------------------------------------

@pytest.fixture(scope="module")
def val_real_alignment():
    """module-scoped: cmalign runs once, result reused across all Val tests.
    skip here rather than via marker so fixture discovery succeeds unconditionally."""
    if not CMALIGN_OK or not CANONICAL_CM:
        pytest.skip("requires cmalign in PATH and SPRINX_CANONICAL_CM")
    seqs = _load_bundle_fa("canonical.fa")
    val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
    aln = sprinx.cmalign_one(val_key, seqs[val_key], CANONICAL_CM)
    return val_key, seqs[val_key], aln


class TestValThreadingFailureIntegration:

    @need_cmalign
    def test_arm_span_check_returns_true_for_real_alignment(self, val_real_alignment):
        """arm_span_has_enough_sequence must return True (threading failure, not
        genuine arm loss) for Val using the real single-seq element span, not the
        multi-seq proxy span from the unit test."""
        _val_key, _val_seq, aln = val_real_alignment
        assert aln is not None, "cmalign_one failed for Val"
        elems = sprinx.get_stem_loop_elements(aln["ss_cons"])
        t_elem = elems[-1]
        assert sprinx.arm_span_has_enough_sequence(aln["aligned_seq"], t_elem), (
            f"arm_span check returned False for Val on real single-seq alignment. "
            f"span={t_elem['span']}, "
            f"n_nts={sum(1 for c in aln['aligned_seq'][t_elem['span'][0]:t_elem['span'][1]] if c not in '-.')}, "
            f"n_stem_cols={len(t_elem['stem_cols'])}"
        )

    @need_cmalign
    def test_patch_produces_balanced_structure_on_real_alignment(self, val_real_alignment):
        """patch_threading_failure_arm must produce a balanced structure when applied
        to the real finalize_structure output with the real T-arm element span.
        this is the exact test that would have caught the production crash: the
        multi-seq proxy in the unit test has different span indices than the real
        single-seq alignment, so the bracket-overlap that caused the crash only
        manifests with the real alignment."""
        _val_key, _val_seq, aln = val_real_alignment
        assert aln is not None
        elems = sprinx.get_stem_loop_elements(aln["ss_cons"])
        t_elem = elems[-1]
        final_seq, final_ss = sprinx.finalize_structure(aln)
        patched = sprinx.patch_threading_failure_arm(
            aln["aligned_seq"], final_seq, final_ss, t_elem
        )
        assert patched.count("(") == patched.count(")"), (
            f"unbalanced patched structure ({patched.count('(')} opens, "
            f"{patched.count(')')} closes):\n{patched}"
        )
        assert len(patched) == len(final_seq), (
            f"patch changed structure length: {len(patched)} vs {len(final_seq)}"
        )

    @need_cmalign
    def test_classify_arm_loss_diagnoses_val_correctly(self, val_real_alignment):
        """classify_arm_loss must diagnose Val as T_OR_VAR_ARM_MISSING (slots 2 and/or 3)
        on the real single-seq alignment, and arm_span_has_enough_sequence must
        return True so it is routed as threading failure, not genuine arm loss."""
        val_key, _val_seq, aln = val_real_alignment
        assert aln is not None
        d = sprinx.classify_arm_loss(val_key, aln["aligned_seq"], aln["ss_cons"])
        assert "T_OR_VAR_ARM_MISSING" in d["call"], (
            f"expected T_OR_VAR_ARM_MISSING call for Val, got {d['call']}"
        )


# -----------------------------------------------------------------------
# select_cm_and_align: end-to-end routing
# -----------------------------------------------------------------------

class TestSelectCmAndAlignIntegration:

    @need_cmalign
    @need_armless
    @pytest.mark.parametrize("fa_key, seq_key", [
        ("canonical_T_human.fa", "mtdbD00063518|Thr|UGU|Homo_sapiens"),
        ("canonical_E_human.fa", "mtdbD00063517|Glu|UUC|Homo_sapiens"),
    ])
    def test_canonical_sequence_not_rerouted(self, fa_key, seq_key):
        seqs = _load_bundle_fa(fa_key)
        armless_index = sprinx.index_armless_cms(ARMLESS_CM_DIR)
        routing = sprinx.select_cm_and_align(seq_key, seqs[seq_key], CANONICAL_CM, armless_index)
        assert routing["final_alignment"] is not None
        assert routing["rerouted"] is False, (
            f"{seq_key}: canonical sequence incorrectly rerouted to armless CM"
        )

    @need_cmalign
    @need_armless
    def test_d_armless_sequence_is_rerouted(self):
        """a genuinely D-armless sequence (human mt-Ser1, GCU anticodon) must be
        rerouted to an armless CM. this test validates the full select_cm_and_align
        pipeline end-to-end with real infernal and real armless CM files."""
        seqs = _load_bundle_fa("D_armless_human.fa")
        seq_key = "mtdbD00063515|Ser1|GCU|Homo_sapiens"
        armless_index = sprinx.index_armless_cms(ARMLESS_CM_DIR)
        routing = sprinx.select_cm_and_align(seq_key, seqs[seq_key], CANONICAL_CM, armless_index)
        assert routing["final_alignment"] is not None
        assert routing["rerouted"] is True, (
            f"D-armless Ser1 was not rerouted (call={routing['diagnosis']['call']})"
        )
        cm_base = os.path.basename(routing["cm_used"])
        assert "wo_d" in cm_base, (
            f"expected a D-armless CM (wo_d), got {cm_base}"
        )

    @need_cmalign
    @need_armless
    def test_val_threading_failure_not_rerouted(self):
        """Val must not be rerouted to an armless CM despite T-arm n_pairs==0:
        arm_span_has_enough_sequence must detect the threading failure and keep
        the canonical CM alignment. this is the end-to-end version of the
        threading failure tests above."""
        seqs = _load_bundle_fa("canonical.fa")
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        armless_index = sprinx.index_armless_cms(ARMLESS_CM_DIR)
        routing = sprinx.select_cm_and_align(val_key, seqs[val_key], CANONICAL_CM, armless_index)
        assert routing["rerouted"] is False, (
"Val incorrectly rerouted to armless CM despite threading failure detection"
        )
        assert routing["threading_failure_elem"] is not None, (
"Val threading failure element not set; arm span check may have failed"
        )


    @need_cmalign
    @need_armless
    def test_doubly_armless_sequence_routes_to_d_and_t_cm(self):
        """a doubly-armless sequence (R. culicivorax mt-Ile, both D and T arms absent;
        Masta & Boore 2008, doi:10.1093/molbev/msn051) must route to a d_and_t CM
        if one exists. classifies as BOTH_ARMS_MISSING (offset=0, D-arm slot and T-arm
        slot both n_pairs==0) rather than UPSTREAM_ARM_MISSING, because with both arms
        absent there is no single-arm register shift."""
        seqs = _load_bundle_fa("both_armless_mature.fa")
        bundle_key = next(k for k in seqs if "culicivorax" in k or "Romanomermis" in k)
        seq = seqs[bundle_key]
        # bundle FASTA uses non-standard field order (id|taxon|aa|anticodon);
        # reformat to id|aa|anticodon|taxon so aa_field_to_cm_code resolves correctly.
        header = "NC_008640.1:3203-3266|Ile|GAU|Romanomermis_culicivorax"
        armless_index = sprinx.index_armless_cms(ARMLESS_CM_DIR)
        d_and_t_available = any(arm == "d_and_t" for _, arm in armless_index)
        routing = sprinx.select_cm_and_align(header, seq, CANONICAL_CM, armless_index)
        diag = routing["diagnosis"]
        assert "BOTH_ARMS_MISSING" in diag["call"] or diag["missing_arm"] in ("d_and_t", "ambiguous"), (
            f"doubly-armless expected BOTH_ARMS_MISSING or d_and_t; got {diag['call']}"
        )
        if d_and_t_available:
            assert routing["rerouted"] is True
            assert "d_and_t" in os.path.basename(routing["cm_used"]), (
                f"expected d_and_t CM; got {routing['cm_used']}"
            )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
