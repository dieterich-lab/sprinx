"""
test_sprinx_unit.py: unit tests for sprinx.label.

no subprocess calls; infernal not required. pre-computed cmalign Stockholm
alignments are in test_data_bundle.txt (==> name <== markers, produced with
`cmalign --notrunc --nonbanded -g`). ViennaRNA (RNA module) is always present
when sprinx is importable (hard import), so RNA API tests are unconditional.

coverage: header parsing, classify_arm_loss, finalize_structure, sprinzl_map
(incl. bulge/junction labeling and the no-unlabeled invariant), arm-span check,
threading patch, CM routing, ViennaRNA API. tests loop over cases internally
rather than parametrizing, to keep the run count small.
"""
import os
import re

import pytest
import RNA

from sprinx import label as sprinx

BUNDLE_PATH = os.path.join(os.path.dirname(__file__), "test_data_bundle.txt")


def load_bundle():
    text = open(BUNDLE_PATH).read()
    chunks = re.split(r"^==> (.+?) <==\n", text, flags=re.MULTILINE)[1:]
    return {name: content for name, content in zip(chunks[0::2], chunks[1::2])}


BUNDLE = load_bundle()


def load_sto(name):
    return sprinx.parse_multi_sto(BUNDLE[name], from_text=True)


def load_fa(name):
    seqs, cur = {}, None
    for line in BUNDLE[name].splitlines():
        if line.startswith(">"):
            cur = line[1:].strip()
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line.strip().upper().replace("T", "U")
    return seqs


def _key(lbl):
    """sort key ordering '23' < '23A' < '23B' < '24', and placing variable-arm
    'e' labels (class-ii: Leu, Ser) between 45 and 46 in true 5'->3' order:
    e11-e17 (5' stem), then e1-e5 (loop), then e21-e27 (3' stem); the 3'
    stem is numbered in reverse as written (e27 comes before e21 reading
    5'->3'), so its rank needs negating or e27 would sort after e21."""
    m = re.match(r"e(\d)(\d)?([a-zA-Z]*)$", lbl)
    if m:
        d1, d2, suffix = m.group(1), m.group(2), m.group(3)
        if d2 is None:
            rank, order = 1, int(d1)
        elif d1 == "1":
            rank, order = 0, int(d2)
        else:
            rank, order = 2, -int(d2)
        return (45, 1, rank, order, suffix)
    m = re.match(r"(\d+)([a-zA-Z]*)", lbl)
    return (int(m.group(1)), 0, 0, 0, m.group(2))


def _monotonic(sprinzl):
    keys = [_key(sprinzl[i]) for i in sorted(sprinzl)]
    return keys == sorted(keys)


def _region(label):
    base = re.match(r"e?\d+", label).group() if label else None
    return sprinx.SPRINZL_REGION.get(base)


def _seq_with_anticodon(ss, stem_index, anticodon):
    """all-A sequence with `anticodon` placed in inner-stem `stem_index`'s loop."""
    topo = sprinx.parse_topology(ss)
    seq = list("A" * len(ss))
    for k, base in enumerate(anticodon):
        seq[topo["inner_stems"][stem_index]["loop_cols"][k]] = base
    return "".join(seq)


# ViennaRNA library contracts (db_from_WUSS, mfe, ptable) are exercised by nearly
# every test below, so they aren't asserted separately. drop_orphan_brackets is a
# defensive net that finalize's both-sides-nulling keeps real data away from, so it
# gets this one direct check on a synthetic orphan.
def test_drop_orphan_brackets_makes_ptable_safe():
    orphaned = "(((((((..((((......)))).(((((.......)))))....................)))))))"
    mangled = orphaned[:5] + ")" + orphaned[6:]   # one unmatched ')'
    assert RNA.ptable(sprinx.drop_orphan_brackets(mangled)) is not None


# -----------------------------------------------------------------------
# header parsing
# -----------------------------------------------------------------------

def test_header_field_extraction():
    cases = [
        ("mtdbD00063517|Glu|UUC|Homo_sapiens", "UUC", "Glu"),
        ("id|Phe|GAT|taxon", "GAU", "Phe"),                       # DNA -> RNA
        ("free text anticodon=GCU here", "GCU", None),            # tag fallback
        # field 3 not a 3nt codon -> None, not silently used as anticodon
        ("NC_008640.1:3203-3266|Romanomermis_culicivorax|Ile|GAU", None, "Romanomermis_culicivorax"),
        ("just a plain header", None, None),
        # GtRNAdb-style 'tRNA-{AA}-{anticodon}' name, no pipes at all
        ("mt-tRNA-Ala-TGC-1-1", "UGC", "Ala"),
        ("mt-tRNA-Leu-TAG-2-1", "UAG", "Leu"),
        ("mt-tRNA-Met-CAT-1-2", "CAU", "Met"),
    ]
    for header, anticodon, aa in cases:
        assert sprinx.header_to_anticodon(header) == anticodon, header
        assert sprinx.header_to_aa(header) == aa, header


# -----------------------------------------------------------------------
# classify_arm_loss
# -----------------------------------------------------------------------

def test_canonical36_no_false_positives():
    seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
    assert len(seqs) == 36
    t_flagged, d_flagged, n_unanchored = [], [], 0
    for name, aligned in seqs.items():
        d = sprinx.classify_arm_loss(name, aligned, ss)
        if d["missing_arm"] == "t":
            t_flagged.append(name)
        if d["missing_arm"] == "d":
            d_flagged.append(name)
        if d["anticodon_search_method"] != "unique_loop_match":
            n_unanchored += 1
    # Val's real T-arm threads into insert columns (n_pairs==0); the arm_span
    # cross-check keeps it canonical. Cys/Ser2 have n_pairs=1 at the D-arm slot
    # under the soft MIN_STEM_PAIRS threshold, likewise resolved by arm_span.
    assert len(t_flagged) <= 1, f"T-arm false positives: {t_flagged}"
    assert len(d_flagged) <= 2, f"D-arm false positives: {d_flagged}"
    assert n_unanchored <= 6, f"{n_unanchored}/36 unanchored"


def test_arm_loss_diagnosis_per_structural_class():
    # canonical: no arm loss, no register shift
    for sto, tag in [("aln_T_canonical_qutrna.sto", "Thr|UGU|Homo"),
                     ("aln_E_canonical_qutrna.sto", "Glu|UUC|Homo"),
                     ("aln_L1_canonical_qutrna.sto", "Leu1|UAG|Homo")]:
        seqs, ss = load_sto(sto)
        name = next(k for k in seqs if tag in k)
        d = sprinx.classify_arm_loss(name, seqs[name], ss)
        assert d["missing_arm"] is None and d["register_offset"] in (0, None), (tag, d["call"])
    # D-armless: register shift offset==1
    seqs, ss = load_sto("aln_S1_qutrna.sto")
    for name, aligned in seqs.items():
        d = sprinx.classify_arm_loss(name, aligned, ss)
        if d["anticodon_stem_index"] is not None:
            assert d["register_offset"] == 1, (name, d["call"])
    # T-armless: never a D-arm call, and no register shift (downstream loss)
    seqs, ss = load_sto("aln_Tarmless_qutrna.sto")
    assert len(seqs) == 17
    for name, aligned in seqs.items():
        d = sprinx.classify_arm_loss(name, aligned, ss)
        assert d["missing_arm"] != "d", (name, d["call"])
        if d["anticodon_stem_index"] is not None:
            assert d["register_offset"] == 0, (name, d["call"])
    # doubly-armless
    seqs, ss = load_sto("aln_both_armless_mature.sto")
    name = "NC_008640.1:3203-3266|Romanomermis_culicivorax|Ile|GAU"
    assert sprinx.classify_arm_loss(name, seqs[name], ss)["missing_arm"] in ("d", "ambiguous", "d_and_t")


def test_d_arm_absent_without_register_shift():
    """no-shift D-arm path (for CMs modeling an extra stem, e.g. TRNAinf-bact.cm).
    no bundled fixture shows it naturally, so blank a real D-arm's own stem columns
    while leaving the anticodon anchor untouched."""
    seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
    name = next(k for k in seqs if "Thr|UGU|Homo" in k)
    base = sprinx.classify_arm_loss(name, seqs[name], ss)
    assert base["missing_arm"] is None, "needs a canonical baseline"
    d_elem = sprinx.get_stem_loop_elements(ss)[base["anticodon_stem_index"] - 1]
    mutated = list(seqs[name])
    for c in d_elem["stem_cols"]:
        mutated[c] = "-"
    mutated = "".join(mutated)
    d = sprinx.classify_arm_loss(name, mutated, ss)
    assert d["register_offset"] == 0 and d["missing_arm"] == "d"
    assert d["call"].startswith("UPSTREAM_ARM_MISSING_slot=")
    assert not sprinx.arm_span_has_enough_sequence(mutated, d_elem)


def test_stem_complementarity_and_anticodon_search():
    # T-armless slot: n_pairs==0 (no column has both partners non-gap)
    seqs, ss = load_sto("aln_Tarmless_qutrna.sto")
    t_elem = sprinx.get_stem_loop_elements(ss)[-1]
    for name, aligned in seqs.items():
        assert sprinx.stem_complementarity(aligned, ss, t_elem)["n_pairs"] == 0, name
    # '.' (multi-seq insert gap) must be stripped alongside '-' before the search,
    # or a '.' in a loop can spuriously match the anticodon.
    seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
    elements = sprinx.get_stem_loop_elements(ss)
    glu = next(k for k in seqs if "Glu|UUC|Homo" in k)
    assert "." in seqs[glu]
    _, method = sprinx.find_anticodon_stem_index(seqs[glu], elements, sprinx.header_to_anticodon(glu))
    assert method == "unique_loop_match"


# -----------------------------------------------------------------------
# sprinzl_map on real data: no unlabeled, monotonic, anticodon placed
# -----------------------------------------------------------------------

def test_sprinzl_map_real_data_invariants():
    # (sto, tag, anticodon, missing_arm)
    cases = [
        ("aln_E_canonical_qutrna.sto",  "Glu|UUC|Homo",  "UUC", None),
        ("aln_T_canonical_qutrna.sto",  "Thr|UGU|Homo",  "UGU", None),
        ("aln_L1_canonical_qutrna.sto", "Leu1|UAG|Homo", "UAG", None),
        # this Ser1 case has an 8-9nt anticodon loop (a stem-edge nucleotide
        # the CM mis-threaded into the loop), not the canonical 7nt: exactly
        # the case _assign_anticodon_loop exists for: the anticodon still
        # lands at 34-35-36 because it's anchored on the anticodon's own
        # position, not the loop's 5' edge.
        ("aln_S1_qutrna.sto",           "Ser1|GCU|Homo", "GCU", "d"),
        ("aln_both_armless_mature.sto", "culicivorax",   "GAU", "d_and_t"),
    ]
    for sto, tag, anticodon, missing_arm in cases:
        seqs, ss_cons = load_sto(sto)
        name = next(k for k in seqs if tag in k)
        seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
        sprinzl = sprinx.sprinzl_map(ss, seq, anticodon, missing_arm)
        assert sprinzl[0] == "1", tag
        assert [i for i in range(len(seq)) if i not in sprinzl] == [], f"{tag}: unlabeled"
        assert _monotonic(sprinzl), f"{tag}: non-monotonic"
        got = "".join(seq[i] for i in sorted(sprinzl) if sprinzl[i] in ("34", "35", "36"))
        assert got == anticodon, f"{tag}: anticodon at 34-36 was {got!r}"
        if missing_arm == "d":
            # option A (Ozerova et al. 2024): the D-armless replacement loop maps
            # onto Sprinzl 8-26 by structural analogy, even with no D-stem pairs.
            labeled = {int(re.match(r"\d+", v).group()) for v in sprinzl.values() if not v.startswith("e")}
            assert any(10 <= p <= 25 for p in labeled), f"{tag}: no D-arm labels"


def test_c_stem3_labels_not_overwritten_by_var_loop():
    """the c-stem-3 columns must keep their own labels (40-43), not get
    overwritten by var_loop's numbering. var_loop's start boundary must be
    the outermost (last) c-stem-3 column, not the innermost one adjacent to
    the loop, or its own assign_slots call bleeds into c-stem-3's columns.
    checked against known-correct values directly, since monotonicity alone
    can't catch a shifted-but-still-increasing run."""
    seqs, ss_cons = load_sto("aln_E_canonical_qutrna.sto")
    name = next(k for k in seqs if "Glu|UUC|Homo" in k)
    seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
    topo = sprinx.parse_topology(ss)
    c_stem3 = topo["inner_stems"][1]["stem3_cols"]
    sprinzl = sprinx.sprinzl_map(ss, seq, "UUC", None)
    assert [sprinzl[p] for p in c_stem3] == ["39", "40", "41", "42", "43"]


def test_variable_arm_stem_gets_e_series_labels():
    """real class-ii case (S. cerevisiae mt-Tyr): a real nested variable-
    arm stem-loop between the c-arm and t-arm gets the Sprinzl e-series
    reserved for it (e11-e17 stem/e1-e5 loop/e21-e27 stem, paired e1N<->e2N),
    and every e-labelled base is complementary to its declared pairing
    partner, not just present."""
    seqs, ss_cons = load_sto("aln_canonical36_qutrna_flags.sto")
    name = next(k for k in seqs if "Tyr|GUA|Saccharomyces" in k)
    seq, ss = sprinx.finalize_structure({"aligned_seq": seqs[name], "ss_cons": ss_cons})
    sprinzl = sprinx.sprinzl_map(ss, seq, "GUA", None)
    assert [i for i in range(len(seq)) if i not in sprinzl] == []
    assert _monotonic(sprinzl)

    by_label = {v: seq[k] for k, v in sprinzl.items()}
    assert {"e11", "e12", "e13", "e1", "e2", "e3", "e4", "e23", "e22", "e21"} <= set(by_label)
    for n in (1, 2, 3):
        pair = (by_label[f"e1{n}"], by_label[f"e2{n}"])
        assert pair in sprinx.WC_PAIRS, f"e1{n}/e2{n} not complementary: {pair}"


def test_assign_anticodon_loop_anchors_on_anticodon_not_loop_edge():
    """direct unit coverage for _assign_anticodon_loop: a loop longer than the
    canonical 7nt (e.g. extra nt on the 5' side) still centres 34-35-36 on
    the real anticodon, with the overflow letter-suffixed onto 33/38 and the
    anticodon's own labels staying fixed."""
    seq = "AAUUGCAUUAA"   # 11nt: 4 before the anticodon, GCA, 4 after
    c_loop = list(range(len(seq)))
    labels = {}
    sprinx._assign_anticodon_loop(labels, seq, c_loop, "GCA")
    assert labels[4] == "34" and labels[5] == "35" and labels[6] == "36"
    assert labels[0] == "32" and labels[1] == "33"
    assert labels[2] == "33A" and labels[3] == "33B"   # overflow before the anchor
    assert labels[7] == "37" and labels[8] == "38"
    assert labels[9] == "38A" and labels[10] == "38B"  # overflow after


def test_long_unpaired_variable_region_uses_e_loop_not_48_overflow():
    """a variable region long enough to signal a real class-ii extension,
    even when no stem was threaded there (so no v_stem is found, e.g. when
    a CM leaves it as one long unpaired run), gets 44/45, e1-e5 (with
    letter-suffix overflow on e5 past 5nt), 46/47/48: the same e-series
    layout as the stem-bearing case, keyed only on length."""
    ss = ("(((((((" + "(((CCC)))" + "(((GGG)))" + "." * 11 + "(((AAA)))" + ")))))))")
    seq = ("A" * 7 + "AAA" + "CCC" + "AAA" + "AAA" + "GGG" + "AAA"
           + "U" * 11 + "AAA" + "AAA" + "AAA" + "A" * 7)
    assert len(ss) == len(seq)
    sprinzl = sprinx.sprinzl_map(ss, seq, "GGG")
    assert [i for i in range(len(seq)) if i not in sprinzl] == []
    var_loop_start = ss.index("." * 11)
    got = [sprinzl[var_loop_start + i] for i in range(11)]
    assert got == ["44", "45", "e1", "e2", "e3", "e4", "e5", "e5A", "46", "47", "48"]


# -----------------------------------------------------------------------
# _forgi_stem_groups: bulge merging vs real arm junctions
# -----------------------------------------------------------------------

def test_forgi_stem_groups():
    classes = {
        "canonical": "(((((((..((((......)))).(((((.......))))).....((((......)))).))))))).",
        "d_armless": "((.((((.......(((((.......)))))....(((((......))))))))))).",
        "t_armless": "(((((((..((((......)))).(((((.......)))))......))))))).",
        "doubly":    "((((............(((((.......))))).........))))....",
    }
    for name, ss in classes.items():
        groups = sprinx._forgi_stem_groups(ss)
        assert sum(1 for g in groups if not g["loop_cols"]) == 1, f"{name}: acceptor count"
    # doubly-armless: acceptor + C-stem joined by an interior loop must not merge
    # (merging leaves no acceptor and makes parse_topology raise).
    assert len(sprinx._forgi_stem_groups(classes["doubly"])) == 2
    topo = sprinx.parse_topology(classes["doubly"])
    assert len(topo["acceptor_5"]) == 4 and len(topo["inner_stems"]) == 1
    # 2nt bulge on each acceptor strand must merge into one acceptor group.
    bulged = "(((..((..((((........)))).(((((.......))))).....(((((.......)))))))..)))."
    topo = sprinx.parse_topology(bulged)
    assert topo["acceptor_5"] == [0, 1, 2, 5, 6]
    assert topo["acceptor_3"] == [65, 66, 69, 70, 71]


# -----------------------------------------------------------------------
# sprinzl_map: bulges/junctions labeled correctly, never unlabeled, monotonic
# -----------------------------------------------------------------------

def test_structural_bulge_labeling():
    acceptor = "(((..((..((((........)))).(((((.......))))).....(((((.......)))))))..)))."
    d_bulge  = "(((((((..(((....))..).(((UUU))).(((....))))))))))"
    t_short  = "(((((((..((((......)))).(((UUU))).....((((......)))).)))))))"

    # (ss, expected {col: region})
    region_checks = [
        # acceptor 5' bulge (3,4) and 3' bulge (67,68) stay acceptor insertions,
        # not absorbed into linker_5 (D-connector) or var_loop.
        (acceptor, {3: "acceptor_5", 4: "acceptor_5", 67: "acceptor_3", 68: "acceptor_3"}),
        # D-stem 3' bulge (18,19) stays a D-stem insertion, not the connector 26.
        (d_bulge, {18: "D_stem_3", 19: "D_stem_3"}),
    ]
    for ss, checks in region_checks:
        sprinzl = sprinx.sprinzl_map(ss, _seq_with_anticodon(ss, 1, "UUU"), "UUU")
        for col, region in checks.items():
            assert _region(sprinzl.get(col, "")) == region, f"col {col} -> {sprinzl.get(col)}"

    # RNAfold-patched 4bp T-stem leaves one unpaired nt (col 52) before the
    # acceptor; it fills the open outermost T-stem-3' slot (65), never blank.
    sprinzl = sprinx.sprinzl_map(t_short, _seq_with_anticodon(t_short, 1, "UUU"), "UUU")
    assert sprinzl.get(52) == "65"

    # every structure: no position left unlabeled, and labels stay monotonic.
    for ss in (acceptor, d_bulge, t_short):
        seq = _seq_with_anticodon(ss, 1, "UUU")
        sprinzl = sprinx.sprinzl_map(ss, seq, "UUU")
        assert [i for i in range(len(seq)) if i not in sprinzl] == [], ss
        assert _monotonic(sprinzl), ss


def test_fill_stem_bulges_overflow_and_ownership():
    # 27 consecutive owned gaps: 27th falls back to a 2-char code, not '['.
    ss = "(" + "." * 27 + ")"
    labels = {0: "5"}
    sprinx._fill_stem_bulges(labels, ss, strands=[[0, len(ss) - 1]])
    assert labels[1] == "5A" and labels[26] == "5Z" and labels[27] == "5AA"
    # a gap outside every strand is not a bulge -> left unlabeled so a bug surfaces.
    labels2 = {1: "10", 6: "20"}
    sprinx._fill_stem_bulges(labels2, ".(....)", strands=[[1, 6]])
    assert 0 not in labels2 and labels2[2] == "10A"


def test_three_stem_double_match_picks_middle_stem():
    # anticodon coincidentally in the D-loop and C-loop of a 3-stem structure:
    # C-stem is the middle by position, and a carried-over missing_arm=t must NOT
    # fire the 2-stem shortcut here.
    ss = "(((((((.(((GGG))).(((GGG))).(((AAA))).)))))))"
    topo = sprinx.parse_topology(ss)
    seq = list("A" * len(ss))
    for stem in (0, 1):
        for k, base in enumerate("GGG"):
            seq[topo["inner_stems"][stem]["loop_cols"][k]] = base
    seq = "".join(seq)
    middle = topo["inner_stems"][1]["stem5_cols"]
    assert sprinx.locate_anticodon_stem(topo, ss, seq, "GGG")["c_stem5"] == middle
    assert sprinx.locate_anticodon_stem(topo, ss, seq, "GGG", "t")["c_stem5"] == middle


# -----------------------------------------------------------------------
# finalize_structure + arm-span/threading-patch
# -----------------------------------------------------------------------

def test_finalize_structure_clean_balanced_equal_length():
    for sto in ("aln_E_canonical_qutrna.sto", "aln_S1_qutrna.sto", "aln_Tarmless_qutrna.sto"):
        seqs, ss_cons = load_sto(sto)
        for name, aligned in seqs.items():
            seq, ss = sprinx.finalize_structure({"aligned_seq": aligned, "ss_cons": ss_cons})
            assert "." not in seq and "-" not in seq, f"{sto}/{name}"
            assert len(seq) == len(ss) and ss.count("(") == ss.count(")"), f"{sto}/{name}"


class TestArmSpanAndPatch:
    """multi-seq bundle proxies the gapped aligned_seq; integration tests exercise
    the same logic on real single-seq cmalign output."""

    CANONICAL_SS = "(((((((..((((......)))).(((((.......)))))....................)))))))"

    def _val(self):
        seqs, ss = load_sto("aln_canonical36_qutrna_flags.sto")
        return seqs, ss, sprinx.get_stem_loop_elements(ss)

    def _synthetic_val_aln(self, seqs, elems):
        base = seqs[next(k for k in seqs if "Val|UAC|Homo" in k)]
        s, e = elems[-1]["span"]
        insert = ("uucaacuuaacuugac" + "-" * (e - s))[:e - s]
        return base[:s] + insert + base[e:]

    def test_span_and_full_span_fold_separate_threading_from_loss(self):
        seqs, ss, elems = self._val()
        t_elem = elems[-1]
        val = next(k for k in seqs if "Val|UAC|Homo" in k)
        # Val's real T-arm passes the span check and folds as a hairpin.
        assert sprinx.arm_span_has_enough_sequence(seqs[val], t_elem)
        assert all(sprinx.arm_span_has_enough_sequence(a, t_elem) for a in seqs.values())
        fseq, _ = sprinx.finalize_structure({"aligned_seq": seqs[val], "ss_cons": ss})
        assert sprinx.arm_is_threading_failure(seqs[val], fseq, t_elem)
        # truly T-armless: fails the span check and does not fold.
        tseqs, tss = load_sto("aln_Tarmless_qutrna.sto")
        te = sprinx.get_stem_loop_elements(tss)[-1]
        assert not any(sprinx.arm_span_has_enough_sequence(a, te) for a in tseqs.values())
        for name, a in tseqs.items():
            fs, _ = sprinx.finalize_structure({"aligned_seq": a, "ss_cons": tss})
            assert not sprinx.arm_is_threading_failure(a, fs, te), name

    def test_patch_recovers_pairs_balanced(self):
        seqs, _, elems = self._val()
        val_seq = load_fa("canonical.fa")[next(k for k in load_fa("canonical.fa") if "Val|UAC|Homo" in k)]
        aln = self._synthetic_val_aln(seqs, elems)
        patched = sprinx.patch_threading_failure_arm("probe", aln, val_seq, self.CANONICAL_SS, elems[-1])
        assert patched.count("(") > self.CANONICAL_SS.count("(")
        assert patched.count("(") == patched.count(")")

    def test_patch_overrides_a_weak_pre_existing_pair_inside_its_own_span(self):
        """a threading-failure span can still carry a single weak pair from
        cmalign's own consensus (below MIN_STEM_PAIRS but still non-'.').
        that pair doesn't block the patch meant to replace it; the patch
        overrides it instead of aborting, since structure inside a confirmed
        threading-failure span is untrustworthy by definition."""
        seqs, _, elems = self._val()
        val_seq = load_fa("canonical.fa")[next(k for k in load_fa("canonical.fa") if "Val|UAC|Homo" in k)]
        aln = self._synthetic_val_aln(seqs, elems)
        # a single bp at the exact position the eventual fold also uses (46-57),
        # directionally consistent with the fold it will be overridden by.
        ss_list = list(self.CANONICAL_SS)
        ss_list[46], ss_list[57] = "(", ")"
        conflicting = "".join(ss_list)
        patched = sprinx.patch_threading_failure_arm("probe", aln, val_seq, conflicting, elems[-1])
        assert patched != conflicting
        assert patched.count("(") == patched.count(")")
        assert patched.count("(") > conflicting.count("(")


# -----------------------------------------------------------------------
# CM library indexing / AA-code routing
# -----------------------------------------------------------------------

def test_cm_filename_indexing(tmp_path):
    for fn in ("armless_trnT_wo_t.cm", "armless_trnL1_wo_d.cm",
               "Metazoa_A.cm", "OtherClade_A.cm", "notacm.txt"):
        (tmp_path / fn).write_text("dummy")
    assert set(sprinx.index_armless_cms(str(tmp_path))) == {("T", "t"), ("L1", "d")}
    canonical = sprinx.index_canonical_cms(str(tmp_path))
    # armless CMs excluded; duplicate aa 'A' keeps exactly one deterministically.
    assert set(canonical) == {"A"}
    assert canonical["A"] in (str(tmp_path / "Metazoa_A.cm"), str(tmp_path / "OtherClade_A.cm"))


def test_aa_field_to_cm_code():
    keys = {("A", "t"), ("V", "d"), ("L1", "t"), ("L2", "d"),
            ("S1", "t"), ("S2", "d"), ("M", "t"), ("T", "t")}
    expected = {"Ala": "A", "Val": "V", "Leu1": "L1", "Leu2": "L2", "Ser1": "S1",
                "Ser2": "S2", "Met": "M", "Xyz": None, "Met3": None,
                # bare aa with no isoacceptor digit (GtRNAdb style): returned
                # as-is when it digit-matches 2+ index entries, for the caller
                # to disambiguate by anticodon instead of failing here.
                "Leu": "L", "Ser": "S"}
    for aa, code in expected.items():
        assert sprinx.aa_field_to_cm_code(aa, keys) == code, aa


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
