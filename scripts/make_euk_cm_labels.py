"""Build tests/data/euk_gtrnadb_cm_labels.tsv from a QutRNA2 seq_to_sprinzl table.

The labels come from QutRNA2's covariance-model path, which predates sprinx and
is what euk users' published coordinates already are. With a QutRNA2 checkout at
$Q, and cmalign on PATH:

    cmalign --notrunc --nonbanded -g -o align.stk \\
        $Q/data/TRNAinf-euk.cm data/cyto/euk_gtrnadb.fa

    python $Q/workflow/scripts/sprinzl_utils.py stk-to-afasta \\
        --output ref.afasta align.stk

    python $Q/workflow/scripts/sprinzl_utils.py consensus-labels \\
        --labels $Q/data/nuclear-euk-masked.txt \\
        --output consensus_labels.tsv align.stk

    python $Q/workflow/scripts/sprinzl_utils.py afasta-to-sprinzl \\
        --consensus-labels consensus_labels.tsv \\
        --output seq_to_sprinzl.tsv ref.afasta

    python scripts/make_euk_cm_labels.py seq_to_sprinzl.tsv \\
        data/cyto/euk_gtrnadb.fa tests/data/euk_gtrnadb_cm_labels.tsv

Those four QutRNA2 commands are what sec_structure.smk cm branch runs, and
the fixture tracks that path.
"""
import collections
import csv
import sys

# the species the QutRNA2 paper reports. Other sequences are dropped, so the
# fixture holds only coordinates somebody published.
SPECIES = ("Homo_sapiens", "Mus_musculus", "Schizosaccharomyces_pombe")

HEADER = """\
# Sprinzl labels for the human, mouse and pombe sequences in
# data/cyto/euk_gtrnadb.fa. QutRNA2's covariance-model path assigned them.
# That path predates sprinx. Moving off it shifts coordinates that users hold.
#
# Written by scripts/make_euk_cm_labels.py - see its docstring to regenerate.
#
# The sequences carry far fewer distinct labelings than there are sequences, so
# each one is written once on a "labels" line, followed by the "seq" lines
# sharing it. Labels run 5'->3', one per base, '-' where the CM assigned none.
"""


def read_fasta(path):
    seqs, name = {}, None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            else:
                seqs[name] += line.strip()
    return seqs


def main(table_path, fasta_path, out_path):
    seqs = {name: seq for name, seq in read_fasta(fasta_path).items()
            if name.startswith(SPECIES)}

    labelled = {}
    with open(table_path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["sprinzl"] == "-" or row["seq_pos"] == "-":
                continue
            labelled[(row["id"].split()[0], int(row["seq_pos"]))] = row["sprinzl"]

    missing = set(seqs) - {name for name, _ in labelled}
    if missing:
        raise SystemExit(f"no labels for {len(missing)} sequences, e.g. {sorted(missing)[:3]}")

    by_pattern = collections.defaultdict(list)
    for name, seq in seqs.items():
        packed = " ".join(labelled.get((name, i), "-") for i in range(1, len(seq) + 1))
        by_pattern[packed].append(name)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(HEADER)
        # commonest labeling first, so the ordinary cases read before the oddities
        for packed, names in sorted(by_pattern.items(), key=lambda kv: (-len(kv[1]), kv[0])):
            fh.write(f"labels\t{packed}\n")
            for name in sorted(names):
                fh.write(f"seq\t{name}\n")

    print(f"{len(by_pattern)} distinct labelings over {len(seqs)} sequences -> {out_path}")


if __name__ == "__main__":
    main(*sys.argv[1:4])