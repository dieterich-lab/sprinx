#!/usr/bin/env python3
"""
convert_output_to_qutrna2-seq_to_sprinzl.py -- convert sprinx.py's sprinzl_mapping.tsv
into QuTRNA2's seq_to_sprinzl.tsv format: one row per (Sprinzl label, tRNA id), giving
that tRNA's 1-indexed sequence position for the label, or '-' if the label doesn't
occur in that particular tRNA.

the reference label set is the union of every distinct sprinzl_position seen across
the whole input, in Sprinzl order -- sprinx.py assigns labels per-sequence (armless
replacement loops, RNAfold-patch overflow, insertion codes), so no fixed master list
exists ahead of time; it has to be built from whatever the input actually contains.

usage: python convert_output_to_qutrna2-seq_to_sprinzl.py sprinzl_mapping.tsv
output: sprinzl_mapping.seq_to_sprinzl.tsv (or --out)
"""
import argparse
import csv
import os
import re
from collections import defaultdict


def _label_sort_key(label):
    """sort key placing every label in true 5'->3' Sprinzl order.
    plain numeric labels (with optional letter-suffix overflow, e.g. '60A')
    sort by (number, suffix length, suffix) -- suffix length before suffix
    itself so a single-letter overflow (Z) sorts before the two-letter
    overflow that follows it (AA), which plain string comparison gets wrong.
    variable-arm 'e' labels (class-ii: Leu, Ser) sit strictly between 45 and
    46: e11-e17 (5' stem) first, then e1-e5 (loop), then e21-e27 (3' stem)
    -- but the 3' stem is numbered in REVERSE as you read it 5'->3' (e27
    comes before e21), so its sort order needs the digit negated or e27
    would wrongly sort after e21."""
    m = re.match(r"e(\d)(\d)?([A-Za-z]*)$", label)
    if m:
        d1, d2, suffix = m.group(1), m.group(2), m.group(3)
        if d2 is None:
            rank, order = 1, int(d1)      # e1..e5: loop
        elif d1 == "1":
            rank, order = 0, int(d2)      # e11..e17: 5' stem
        else:
            rank, order = 2, -int(d2)     # e21..e27: 3' stem, reverse order
        return (45, 1, rank, order, len(suffix), suffix)
    m = re.match(r"(\d+)([A-Za-z]*)", label)
    num, suffix = m.group(1), m.group(2)
    return (int(num), 0, 0, 0, len(suffix), suffix)


def convert(in_path, out_path):
    by_id = defaultdict(dict)   # seq_id -> {label: seq_pos}
    order = []                  # seq_id in order of first appearance
    labels = set()

    with open(in_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            label = row["sprinzl_position"].strip()
            if not label:            # unlabeled position (see sprinx.py's own
                continue              # unlabeled-position warning); nothing to map
            seq_id = row["seq_id"]
            if seq_id not in by_id:
                order.append(seq_id)
            # qutrna2 convention: 17A/20A/20B, not 17a/20a/20b -- but variable-arm
            # 'e' labels (e11, e21, ...) keep their lowercase 'e', the standard
            # convention for them, untouched by this historical-code uppercasing.
            if not label.startswith(("e", "E")):
                label = label.upper()
            by_id[seq_id][label] = int(row["seq_index"]) + 1   # 1-indexed seq_pos
            labels.add(label)

    ordered_labels = sorted(labels, key=_label_sort_key)

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["sprinzl", "seq_pos", "id"])
        for seq_id in order:
            w.writerow(["-1", "-", seq_id])
            positions = by_id[seq_id]
            for label in ordered_labels:
                w.writerow([label, positions.get(label, "-"), seq_id])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sprinzl_tsv", help="sprinx.py output TSV (sprinzl_mapping.tsv)")
    parser.add_argument("--out", default=None,
                        help="output path (default: <input>.seq_to_sprinzl.tsv)")
    args = parser.parse_args()

    out_path = args.out or os.path.splitext(args.sprinzl_tsv)[0] + ".seq_to_sprinzl.tsv"
    convert(args.sprinzl_tsv, out_path)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
