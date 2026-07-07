"""
test_sprinx_unit.py -- unit tests for sprinx.py.

no subprocess calls; infernal not required. pre-computed cmalign Stockholm
alignments are in test_data_bundle.txt (==> name <== markers, produced with
`cmalign --notrunc --nonbanded -g`). ViennaRNA (RNA module) is always present
when sprinx is importable (hard import), so RNA API tests are unconditional.

coverage:
  header parsing, classify_arm_loss, finalize_structure, sprinzl_map,
  arm-span check, patch pre/post-validation, CM routing, ViennaRNA API.
"""
import os
import re
import sys

import pytest
import RNA

sys.path.insert(0, os.path.dirname(__file__))
import sprinx

BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "test_data_bundle.txt")


def load_bundle():
    text = open(BUNDLE_PATH).read()
    chunks = re.split(r"^==> (.+?) <==\n", text, flags=re.MULTILINE)[1:]
    return {name: content for name, content in zip(chunks[0::2], chunks[1::2])}


BUNDLE = load_bundle()


def load_sto(name):
    return sprinx.parse_multi_sto(BUNDLE[name], from_text=True)


def load_fa(name):
    """parse a FASTA block from the bundle into {header: seq}."""
    seqs = {}
    cur = None
    for line in BUNDLE[name].splitlines():
        if line.startswith(">"):
            cur = line[1:].strip()
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip().upper().replace("T", "U")
    return seqs


# -----------------------------------------------------------------------
# ViennaRNA API (unconditional: RNA always importable alongside sprinx)
# -----------------------------------------------------------------------

class TestViennaRNA:

    def test_db_from_wuss_same_length(self):
        """RNA.db_from_WUSS must return same-length string; only ()."""
        wuss = "(((((((,,<<<<______>>>>,<<<<<_______>>>>>,,,,..,,,...........))))))):"
        db = RNA.db_from_WUSS(wuss)
        assert len(db) == len(wuss)
        assert all(c in "()." for c in db)

    def test_db_from_wuss_idempotent_on_dotbracket(self):
        """db_from_WUSS on plain dot-bracket must return it unchanged."""
        db = "(((..(((....))).)))"
        assert RNA.db_from_WUSS(db) == db

    def test_mfe_returns_balanced_structure(self):
        """fold_compound.mfe() on a known hairpin-former must return a
        balanced dot-bracket with at least one pair. sequence is the Val
        T-arm subseq extracted in threading patch tests."""
        ss, mfe = RNA.fold_compound("UUCAACUUAACUUGAC").mfe()
        assert len(ss) == 16
        assert ss.count("(") == ss.count(")")
        assert "(" in ss

    def test_pair_table_raises_on_unbalanced_input(self):
        """pair_table must raise ValueError on unbalanced brackets; this is
        the exception path that patch_threading_failure_arm pre-validation
        was added to prevent reaching."""
        with pytest.raises(ValueError, match="ViennaRNA rejected"):
            sprinx.pair_table("((((((....)")  # 6 opens, 1 close

    def test_drop_orphan_brackets_makes_ptable_safe(self):
        """an orphan bracket from gap-stripping must be made inert before
        RNA.ptable sees it; drop_orphan_brackets + pair_table must not raise."""
        orphaned = "(((((((..((((......)))).(((((.......)))))....................)))))))"
        # introduce one unmatched ')' that would unbalance the structure
        mangled = orphaned[:5] + ")" + orphaned[6:]
        fixed = sprinx.drop_orphan_brackets(mangled)
        pairs = sprinx.pair_table(fixed)
        assert isinstance(pairs, dict)


# -----------------------------------------------------------------------
# header parsing
# -----------------------------------------------------------------------

class TestHeaderParsing:

    def test_standard_pipe_format(self):
        assert sprinx.header_to_anticodon("mtdbD00063517|Glu|UUC|Homo_sapiens") == "UUC"
        assert sprinx.header_to_aa("mtdbD00063517|Glu|UUC|Homo_sapiens") == "Glu"

    def test_dna_alphabet_converted_to_rna(self):
        assert sprinx.header_to_anticodon("id|Phe|GAT|taxon") == "GAU"

    def test_anticodon_equals_tag_fallback(self):
        assert sprinx.header_to_anticodon("free text anticodon=GCU here") == "GCU"

    def test_mismatched_field_order_fails_safe(self):
        # field 3 is not 3nt -> None, not silently used as anticodon
        assert sprinx.header_to_anticodon("NC_008640.1:3203-3266|Romanomermis_culicivorax|Ile|GAU") is None

    def test_no_anticodon_anywhere(self):
        assert sprinx.header_to_anticodon("just a plain header") is None


# -----------------------------------------------------------------------
# canonical-36: Homo, Saccharomyces, Schizosaccharomyces. zero false positives.
# -----------------------------------------------------------------------

class TestCanonical36NoFalsePositives:

    @pytest.fixture(scope="class")
    def alignment(self):
        return load_sto("aln_canonical36_qutrna_flags.sto")

    def test_36_sequences_present(self, alignment):
        seqs, ss = alignment
        assert len(seqs) == 36

    def test_no_false_positives(self, alignment):
        seqs, ss = alignment
        n_unanchored = 0
        t_arm_flagged = []
        d_arm_flagged = []
        for name, aligned_seq in seqs.items():
            d = sprinx.classify_arm_loss(name, aligned_seq, ss)
            if d["missing_arm"] == "d":
                d_arm_flagged.append((name, d["call"]))
            if d["missing_arm"] == "t":
                t_arm_flagged.append((name, d["call"]))
            if d["anticodon_search_method"] != "unique_loop_match":
                n_unanchored += 1
        # human mt-Val has n_pairs==0 at T-arm slot in this multi-seq alignment
        # (cryo-EM confirms real T-arm; Suzuki et al. 2020 PMID:32075901). tolerance 1.
        assert len(t_arm_flagged) <= 1, f"unexpected T-arm false positives: {t_arm_flagged}"
        # D-arm loss detected via register shift (offset>0) alone has zero false
        # positives here (the original, still-true invariant). classify_arm_loss
        # ALSO flags D-arm loss with no shift (offset==0, D-arm's own slot under
        # MIN_STEM_PAIRS) -- needed for CMs that model an extra stem (see
        # UPSTREAM_ARM_MISSING_slot=, added to catch a real bug in TRNAinf-bact.cm
        # routing). that second path uses the same soft MIN_STEM_PAIRS threshold
        # as the T-arm check, so it can occasionally flag a real-but-weakly-paired
        # D-arm too (mtdbD00063507|Cys and mtdbD00063509|Ser2 in this fixture:
        # n_pairs=1 at the D-arm slot, both otherwise fully canonical).
        # classify_arm_loss alone doesn't resolve that ambiguity -- only
        # select_cm_and_align's arm_span_has_enough_sequence cross-check does,
        # which correctly identifies both as "enough sequence, not genuine loss"
        # and keeps them on the canonical CM (verified separately, not rerouted).
        # tolerance 2, same reason as the T-arm tolerance above.
        assert len(d_arm_flagged) <= 2, f"unexpected D-arm false positives: {d_arm_flagged}"
        # small rate of unanchored (AT-rich anticodons) expected; >6 suggests regression
        assert n_unanchored <= 6, f"{n_unanchored}/36 fell to ambiguous-anchor fallback"


# -----------------------------------------------------------------------
# individual canonical (human Thr/Glu/Leu1)
# -----------------------------------------------------------------------

class TestIndividualCanonical:

    @pytest.mark.parametrize("sto_name, seq_name", [
        ("aln_T_canonical_qutrna.sto",  "mtdbD00063518|Thr|UGU|Homo_sapiens"),
        ("aln_E_canonical_qutrna.sto",  "mtdbD00063517|Glu|UUC|Homo_sapiens"),
        ("aln_L1_canonical_qutrna.sto", "mtdbD00063516|Leu1|UAG|Homo_sapiens"),
    ])
    def test_no_arm_loss(self, sto_name, seq_name):
        seqs, ss = load_sto(sto_name)
        d = sprinx.classify_arm_loss(seq_name, seqs[seq_name], ss)
        assert d["missing_arm"] is None, f"{seq_name}: {d['call']}"
        assert d["register_offset"] in (0, None)


# -----------------------------------------------------------------------
# T-armless: 17 seqs (Ascaris suum, Habronattus oregonensis)
# -----------------------------------------------------------------------

class TestTArmless:

    @pytest.fixture(scope="class")
    def alignment(self):
        return load_sto("aln_Tarmless_qutrna.sto")

    def test_17_sequences_present(self, alignment):
        seqs, ss = alignment
        assert len(seqs) == 17

    def test_d_arm_never_flagged_offset_always_zero_when_anchored(self, alignment):
        seqs, ss = alignment
        for name, aligned_seq in seqs.items():
            d = sprinx.classify_arm_loss(name, aligned_seq, ss)
            assert d["missing_arm"] != "d", f"{name}: false D-arm call ({d['call']})"
            if d["anticodon_stem_index"] is not None:
                assert d["register_offset"] == 0, (
                    f"{name}: T-arm loss is downstream; must not shift register "
                    f"(offset={d['register_offset']})"
                )


# -----------------------------------------------------------------------
# D-armless: 3 seqs (Sphenodon, Mus, Homo Ser1)
# -----------------------------------------------------------------------

class TestDArmless:

    @pytest.fixture(scope="class")
    def alignment(self):
        return load_sto("aln_S1_qutrna.sto")

    def test_3_d_armless_sequences_all_have_offset_1(self, alignment):
        seqs, ss = alignment
        for name, aligned_seq in seqs.items():
            d = sprinx.classify_arm_loss(name, aligned_seq, ss)
            if d["anticodon_stem_index"] is not None:
                assert d["register_offset"] == 1, (
                    f"{name}: D-armless expected offset=1, got {d['register_offset']} ({d['call']})"
                )


# -----------------------------------------------------------------------
# D-arm loss with NO register shift: the path added for CMs (e.g. a bacterial
# whole-family CM) that model more than the canonical D/C/T trio, where D-arm
# loss doesn't always produce the offset>0 shift the qutrna fixtures above
# always show. no bundled fixture naturally exercises this (all bundled D-arm
# loss is via register shift), so it's reproduced by blanking a real, present
# D-arm's own stem columns on an otherwise-canonical sequence -- confirmed
# during development this reproduces the real bacterial-CM bug exactly.
# -----------------------------------------------------------------------

class TestDArmNoShiftDetection:

    @pytest.fixture(scope="class")
    def canonical_data(self):
        return load_sto("aln_canonical36_qutrna_flags.sto")

    def _blank_d_arm(self, aligned_seq, ss, diagnosis):
        """gap out a real D-arm's own stem columns, forcing absent() to flag it
        while leaving the anticodon anchor (and thus register_offset) untouched."""
        d_elem = sprinx.get_stem_loop_elements(ss)[diagnosis["anticodon_stem_index"] - 1]
        mutated = list(aligned_seq)
        for c in d_elem["stem_cols"]:
            mutated[c] = "-"
        return "".join(mutated)

    def test_d_arm_absent_no_shift_is_flagged(self, canonical_data):
        seqs, ss = canonical_data
        name = next(k for k in seqs if "Thr|UGU|Homo" in k)
        baseline = sprinx.classify_arm_loss(name, seqs[name], ss)
        assert baseline["missing_arm"] is None, "test needs a genuinely canonical baseline"
        mutated = self._blank_d_arm(seqs[name], ss, baseline)
        d = sprinx.classify_arm_loss(name, mutated, ss)
        assert d["register_offset"] == 0, "blanking stem columns must not move the anchor"
        assert d["missing_arm"] == "d"
        assert d["call"].startswith("UPSTREAM_ARM_MISSING_slot="), d["call"]

    def test_d_arm_absent_no_shift_span_check_says_genuine_loss(self, canonical_data):
        """the blanked D-arm has no sequence left in its span at all, so it must
        fail arm_span_has_enough_sequence -- select_cm_and_align's cross-check
        for this path (see its docstring) would correctly reroute, not patch."""
        seqs, ss = canonical_data
        name = next(k for k in seqs if "Thr|UGU|Homo" in k)
        baseline = sprinx.classify_arm_loss(name, seqs[name], ss)
        mutated = self._blank_d_arm(seqs[name], ss, baseline)
        d = sprinx.classify_arm_loss(name, mutated, ss)
        d_elem = sprinx.get_stem_loop_elements(ss)[d["anticodon_stem_index"] - 1]
        assert not sprinx.arm_span_has_enough_sequence(mutated, d_elem)


# -----------------------------------------------------------------------
# doubly-armless: Romanomermis culicivorax Ile (mature)
# -----------------------------------------------------------------------

class TestBothArmlessMature:

    @pytest.fixture(scope="class")
    def alignment(self):
        return load_sto("aln_both_armless_mature.sto")

    def test_both_arms_structurally_absent(self, alignment):
        seqs, ss = alignment
        name = "NC_008640.1:3203-3266|Romanomermis_culicivorax|Ile|GAU"
        d = sprinx.classify_arm_loss(name, seqs[name], ss)
        assert d["missing_arm"] in ("d", "ambiguous", "d_and_t"), (
            f"expected d, ambiguous, or d_and_t; got {d['missing_arm']} ({d['call']})"
        )


# -----------------------------------------------------------------------
# stem_complementarity
# -----------------------------------------------------------------------

class TestStemComplementarity:

    def test_structurally_impossible_stem_has_zero_pairs(self):
        """T-armless sequences must have n_pairs==0 at the T-arm slot: if no
        alignment column has both pairing partners non-gap, no stem can exist,
        regardless of nucleotide content."""
        seqs, ss = load_sto("aln_Tarmless_qutrna.sto")
        t_elem = sprinx.get_stem_loop_elements(ss)[-1]
        for name, aligned_seq in seqs.items():
            result = sprinx.stem_complementarity(aligned_seq, ss, t_elem)
            assert result["n_pairs"] == 0, (
                f"{name}: T-armless has n_pairs={result['n_pairs']} at T-arm slot"
            )

    def test_two_gap_symbols_both_filtered(self):
        """'.' (multi-seq insert gap) must be stripped from loop sequences alongside
        '-'. without dot filtering, a '.' at a loop position can spuriously match the
        anticodon search, producing an ambiguous result. tested on Glu|UUC in the
        canonical-36 multi-seq alignment (has '.' chars; UUC is rare in loop sequences
        so any spurious match would be immediately visible)."""
        seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
        elements = sprinx.get_stem_loop_elements(ss)
        glu_key = next(k for k in seqs if "Glu|UUC|Homo" in k)
        assert "." in seqs[glu_key], "test requires '.' in aligned_seq (multi-seq context)"
        idx, method = sprinx.find_anticodon_stem_index(
            seqs[glu_key], elements, sprinx.header_to_anticodon(glu_key)
        )
        assert method == "unique_loop_match", (
            f"Glu|UUC should anchor uniquely; got {method!r}. "
            f"ambiguous result suggests '.' not stripped from loop_seq."
        )


# -----------------------------------------------------------------------
# Sprinzl assignment
# -----------------------------------------------------------------------

class TestSprinzlAssignment:

    def test_canonical_seq_sprinzl_pos_1_at_first_nt(self):
        seqs, ss_cons = load_sto("aln_E_canonical_qutrna.sto")
        name = "mtdbD00063517|Glu|UUC|Homo_sapiens"
        seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
        sprinzl = sprinx.sprinzl_map(ss, seq, "UUC")
        assert sprinzl[0] == "1"

    def test_anticodon_lands_in_c_loop(self):
        seqs, ss_cons = load_sto("aln_E_canonical_qutrna.sto")
        name = "mtdbD00063517|Glu|UUC|Homo_sapiens"
        seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
        sprinzl = sprinx.sprinzl_map(ss, seq, "UUC")
        found = any(
            len(regions) == 3 and all(r == "C_loop" for r in regions)
            for ac_pos in [m.start() for m in re.finditer("(?=UUC)", seq)]
            for regions in [[sprinx.SPRINZL_REGION.get(re.match(r"\d+", lbl).group())
                             for lbl in [sprinzl.get(p) for p in range(ac_pos, ac_pos + 3)]
                             if lbl]]
        )
        assert found, "no UUC occurrence landed in C_loop region"

    def test_d_armless_replacement_loop_gets_d_arm_labels(self):
        """option A convention (Ozerova et al. 2024, PMC11571959): D-armless tRNAs
        have a replacement loop in place of the D-arm; map it onto Sprinzl positions
        8-26 by structural analogy. the D-arm stem position slots (10-13, 22-25) are
        assigned to linker nucleotides even though no stem pairs exist there."""
        seqs, ss_cons = load_sto("aln_S1_qutrna.sto")
        name = "mtdbD00063515|Ser1|GCU|Homo_sapiens"
        seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
        sprinzl = sprinx.sprinzl_map(ss, seq, "GCU")
        assigned_regions = {sprinx.SPRINZL_REGION.get(re.match(r"\d+", lbl).group())
                            for lbl in sprinzl.values() if lbl}
        # replacement loop must receive at least one D-arm region label beyond
        # the connector positions (8-9). D_stem_5 is sufficient: not all replacement
        # loops are long enough to reach D_loop (pos 14), and that is biologically valid.
        assert any(r in assigned_regions for r in ("D_stem_5", "D_loop", "D_stem_3")), (
            "D-armless linker must receive D-arm Sprinzl labels beyond connector (option A; PMC11571959)"
        )
        # verify numerically: at least one position in 10-25 is labeled
        d_arm_labeled = {lbl for lbl in sprinzl.values()
                         if lbl and re.match(r"\d+", lbl)
                         and 10 <= int(re.match(r"\d+", lbl).group()) <= 25}
        assert d_arm_labeled, "no positions 10-25 assigned to D-armless replacement loop"


# -----------------------------------------------------------------------
# finalize_structure: gap symbol stripping invariants
# -----------------------------------------------------------------------

class TestFinalizeStructure:

    def test_strips_dot_gap_symbols_from_multi_seq_aligned(self):
        """in multi-seq cmalign output, '.' marks a gap at another sequence's insert
        position -- not real sequence for this record. must be stripped alongside '-'."""
        seqs, ss_cons = load_sto("aln_canonical36_qutrna_flags.sto")
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[val_key], "ss_cons": ss_cons})
        assert "." not in seq
        assert "-" not in seq
        expected_len = sum(1 for c in seqs[val_key] if c not in "-.")
        assert len(seq) == expected_len
        assert len(ss) == len(seq)

    def test_seq_and_ss_same_length_across_alignment_types(self):
        """seq and ss from finalize_structure must always be equal length."""
        for sto_name in ("aln_E_canonical_qutrna.sto", "aln_S1_qutrna.sto",
                         "aln_Tarmless_qutrna.sto"):
            seqs, ss_cons = load_sto(sto_name)
            for name, aligned_seq in seqs.items():
                seq, ss = sprinx.finalize_structure({"aligned_seq": aligned_seq, "ss_cons": ss_cons})
                assert len(seq) == len(ss), (
                    f"{sto_name}/{name}: seq ({len(seq)}) and ss ({len(ss)}) lengths differ"
                )

    def test_output_ss_is_balanced(self):
        """drop_orphan_brackets in finalize_structure must ensure balanced output
        even when gap-stripping removes one side of a paired column."""
        for sto_name in ("aln_E_canonical_qutrna.sto", "aln_S1_qutrna.sto"):
            seqs, ss_cons = load_sto(sto_name)
            for name, aligned_seq in seqs.items():
                _, ss = sprinx.finalize_structure({"aligned_seq": aligned_seq, "ss_cons": ss_cons})
                assert ss.count("(") == ss.count(")"), (
                    f"{sto_name}/{name}: unbalanced ss from finalize_structure"
                )


# -----------------------------------------------------------------------
# arm-span sequence check and threading patch (using multi-seq bundle proxy)
# -----------------------------------------------------------------------

class TestArmSpanAndPatch:
    """the multi-seq bundle is used as a proxy for the gapped aligned_seq here.
    integration tests (test_sprinx_integration.py) verify the same logic with
    single-seq cmalign output against real CMs, which is the production code path."""

    @pytest.fixture(scope="class")
    def canonical36_data(self):
        seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
        return seqs, sprinx.get_stem_loop_elements(ss)

    @pytest.fixture(scope="class")
    def tarmless_data(self):
        seqs, ss = load_sto("aln_Tarmless_qutrna.sto")
        return seqs, sprinx.get_stem_loop_elements(ss)

    def test_val_detected_as_threading_failure(self, canonical36_data):
        """human mt-Val has a real T-arm (cryo-EM, Suzuki et al. 2020 PMID:32075901)
        that goes into insert columns under the canonical CM; span check must
        return True (threading failure, not genuine arm loss)."""
        seqs, elems = canonical36_data
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        assert sprinx.arm_span_has_enough_sequence(seqs[val_key], elems[-1])

    def test_all_canonical_sequences_pass_span_check(self, canonical36_data):
        seqs, elems = canonical36_data
        t_elem = elems[-1]
        fails = [n for n, a in seqs.items() if not sprinx.arm_span_has_enough_sequence(a, t_elem)]
        assert fails == [], f"spurious arm-loss from span check: {fails}"

    def test_all_t_armless_sequences_fail_span_check(self, tarmless_data):
        seqs, elems = tarmless_data
        t_elem = elems[-1]
        passes = [n for n, a in seqs.items() if sprinx.arm_span_has_enough_sequence(a, t_elem)]
        assert passes == [], f"T-armless mis-classified as threading failure: {passes}"

    def test_val_arm_is_threading_failure_via_full_span(self, canonical36_data):
        """regression test: a real arm's sequence isn't always in insert (lowercase)
        columns -- human mt-Val's T-arm sits partly in matched columns here.
        arm_is_threading_failure and patch_threading_failure_arm both use the FULL
        span (matched + insert together, _arm_full_span_subseq_and_fold), so they see
        the same hairpin. an earlier version used an insert-only extraction for
        detection that found nothing here, wrongly causing Val to be rerouted instead
        of patched."""
        seqs, elems = canonical36_data
        _, ss = load_sto("aln_canonical36_qutrna_flags.sto")
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        aligned_seq = seqs[val_key]
        t_elem = elems[-1]
        final_seq, _ = sprinx.finalize_structure({"aligned_seq": aligned_seq, "ss_cons": ss})
        assert sprinx.arm_is_threading_failure(aligned_seq, final_seq, t_elem)

    def test_val_patch_recovers_real_data_end_to_end(self, canonical36_data):
        """regression test: patch_threading_failure_arm must actually pair up real
        (unmodified) human mt-Val data, not just the synthetic aln used by
        test_patch_recovers_paired_structure. an earlier version silently no-opped
        here because it sourced the patch from insert-only columns while Val's T-arm
        sequence sits mostly in matched columns -- detection said "real hairpin" but
        the patch itself found nothing to splice in, so the final structure kept the
        T-arm as unpaired dots despite the "patched" label."""
        seqs, elems = canonical36_data
        _, ss = load_sto("aln_canonical36_qutrna_flags.sto")
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        aligned_seq = seqs[val_key]
        t_elem = elems[-1]
        final_seq, final_ss = sprinx.finalize_structure({"aligned_seq": aligned_seq, "ss_cons": ss})
        patched = sprinx.patch_threading_failure_arm(aligned_seq, final_seq, final_ss, t_elem)
        assert patched.count("(") > final_ss.count("("), "patch did not add paired positions"
        assert patched.count("(") == patched.count(")")

    def test_t_armless_sequences_fail_full_span_check(self, tarmless_data):
        """complement to test_all_t_armless_sequences_fail_span_check: genuinely
        T-armless sequences must not fold as a hairpin under the full-span check either."""
        seqs, elems = tarmless_data
        t_elem = elems[-1]
        _, ss = load_sto("aln_Tarmless_qutrna.sto")
        folds = []
        for name, aligned_seq in seqs.items():
            final_seq, _ = sprinx.finalize_structure({"aligned_seq": aligned_seq, "ss_cons": ss})
            if sprinx.arm_is_threading_failure(aligned_seq, final_seq, t_elem):
                folds.append(name)
        assert folds == [], f"T-armless sequence folded as a real hairpin: {folds}"

    @staticmethod
    def _synthetic_val_aln(seqs, elems):
        """single-seq cmalign puts t-arm inserts as lowercase chars within the t-arm
        element span; multi-seq has them as uppercase model matches or outside the span,
        so c.islower() returns nothing there. replace span content with known t-arm
        sequence (lowercase) + dashes while keeping the prefix char count for ungapped_idx."""
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        base = seqs[val_key]
        s, e = elems[-1]["span"]
        insert = ("uucaacuuaacuugac" + "-" * (e - s))[:e - s]
        return base[:s] + insert + base[e:]

    def test_patch_recovers_paired_structure(self, canonical36_data):
        """c.islower() fix: only insert chars (lowercase) are extracted from the span.
        multi-seq Val has 4 upstream uppercase model-match chars within the span that
        dilute arm_subseq to AGAUUUCAACUUAAC (mfe=0); single-seq-like synthetic aln
        gives UUCAACUUAACUUGAC which folds to .((((......))))."""
        seqs, elems = canonical36_data
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        val_seq = load_fa("canonical.fa")[next(k for k in load_fa("canonical.fa") if "Val|UAC|Homo" in k)]
        aln = self._synthetic_val_aln(seqs, elems)
        canonical_ss = "(((((((..((((......)))).(((((.......)))))....................)))))))"
        patched = sprinx.patch_threading_failure_arm(aln, val_seq, canonical_ss, elems[-1])
        assert patched.count("(") > canonical_ss.count("("), "patch did not add paired positions"

    def test_patch_produces_balanced_structure(self, canonical36_data):
        seqs, elems = canonical36_data
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        val_seq = load_fa("canonical.fa")[next(k for k in load_fa("canonical.fa") if "Val|UAC|Homo" in k)]
        aln = self._synthetic_val_aln(seqs, elems)
        canonical_ss = "(((((((..((((......)))).(((((.......)))))....................)))))))"
        patched = sprinx.patch_threading_failure_arm(aln, val_seq, canonical_ss, elems[-1])
        assert patched.count("(") == patched.count(")"), (
            f"unbalanced patch: {patched.count('(')} opens, {patched.count(')')} closes"
        )

    def test_patch_aborts_on_bracket_overlap(self, canonical36_data):
        """pre-validation aborts when arm_ss brackets would overwrite non-dot positions.
        arm_ss = .((((......)))).  -> opens at ungapped_positions[1..4]=[46..49],
        closes at [56..59]. placing ) at 46-48 triggers abort."""
        seqs, elems = canonical36_data
        val_key = next(k for k in seqs if "Val|UAC|Homo" in k)
        val_seq = load_fa("canonical.fa")[next(k for k in load_fa("canonical.fa") if "Val|UAC|Homo" in k)]
        aln = self._synthetic_val_aln(seqs, elems)
        ss_list = list("(((((((..((((......)))).(((((.......)))))....................)))))))")
        ss_list[46] = ")"
        ss_list[47] = ")"
        ss_list[48] = ")"
        original = "".join(ss_list)
        patched = sprinx.patch_threading_failure_arm(aln, val_seq, original, elems[-1])
        assert patched == original, "patch must abort when arm_ss brackets overlap existing ')'"



# -----------------------------------------------------------------------
# CM library indexing / AA-code routing
# -----------------------------------------------------------------------

class TestCmRouting:

    def test_armless_cm_filename_parsing(self, tmp_path):
        (tmp_path / "armless_trnT_wo_t.cm").write_text("dummy")
        (tmp_path / "armless_trnL1_wo_d.cm").write_text("dummy")
        (tmp_path / "TRNAinf-euk.cm").write_text("dummy")
        index = sprinx.index_armless_cms(str(tmp_path))
        assert set(index.keys()) == {("T", "t"), ("L1", "d")}

    def test_canonical_cm_filename_parsing(self, tmp_path):
        (tmp_path / "Metazoa_A.cm").write_text("dummy")
        (tmp_path / "Metazoa_L1.cm").write_text("dummy")
        (tmp_path / "armless_trnP_wo_d.cm").write_text("dummy")  # must be excluded
        (tmp_path / "notacm.txt").write_text("dummy")
        index = sprinx.index_canonical_cms(str(tmp_path))
        assert set(index.keys()) == {"A", "L1"}

    def test_canonical_cm_duplicate_aa_keeps_one_deterministically(self, tmp_path):
        (tmp_path / "Metazoa_A.cm").write_text("dummy")
        (tmp_path / "OtherClade_A.cm").write_text("dummy")
        index = sprinx.index_canonical_cms(str(tmp_path))
        assert list(index.keys()) == ["A"]
        assert index["A"] in (str(tmp_path / "Metazoa_A.cm"), str(tmp_path / "OtherClade_A.cm"))

    def test_aa_field_to_cm_code_standard(self):
        keys = {("A", "t"), ("V", "d"), ("T", "t")}
        assert sprinx.aa_field_to_cm_code("Ala", keys) == "A"
        assert sprinx.aa_field_to_cm_code("Val", keys) == "V"
        assert sprinx.aa_field_to_cm_code("Xyz", keys) is None

    def test_aa_field_to_cm_code_isoacceptors(self):
        keys = {("L1", "t"), ("L2", "d"), ("S1", "t"), ("S2", "d"), ("M", "t")}
        assert sprinx.aa_field_to_cm_code("Leu1", keys) == "L1"
        assert sprinx.aa_field_to_cm_code("Leu2", keys) == "L2"
        assert sprinx.aa_field_to_cm_code("Ser1", keys) == "S1"
        assert sprinx.aa_field_to_cm_code("Ser2", keys) == "S2"
        assert sprinx.aa_field_to_cm_code("Met",  keys) == "M"
        assert sprinx.aa_field_to_cm_code("Met3", keys) is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))