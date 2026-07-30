#!/usr/bin/env python3
"""
generate_synthetic_cyto_seqs.py: rebuild data/cyto/{euk,arch,bact}.fa from
the bundled isotype CM databases: cmfetch one CM per amino acid, cmemit -c
for its consensus sequence, then derive an anticodon by aligning the
consensus back to its own CM and taking the middle 3nt of the middle
stem-loop (the anticodon loop in a canonical cloverleaf).

usage: python scripts/generate_synthetic_cyto_seqs.py
"""
import os
import tempfile

from sprinx.common import cmalign_one, get_stem_loop_elements, run

DOMAIN_AAS = {
    "euk": ["Ala", "Leu", "SeC"],
    "arch": ["Ala", "Leu"],
    "bact": ["Ala", "Leu"],
}


def consensus_seq(cm_path):
    stdout, stderr, rc = run(["cmemit", "-c", cm_path])
    if rc != 0:
        raise RuntimeError(f"cmemit -c failed: {stderr.strip()}")
    return "".join(stdout.splitlines()[1:]).upper()


def derive_anticodon(seq, cm_path):
    aln = cmalign_one("consensus", seq, cm_path)
    elems = get_stem_loop_elements(aln["ss_cons"])
    loop = "".join(aln["aligned_seq"][c] for c in elems[1]["loop_cols"]
                   if aln["aligned_seq"][c] not in "-.")
    mid = len(loop) // 2
    return loop[mid - 1:mid + 2]


def main():
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    for domain, aas in DOMAIN_AAS.items():
        cm_db = f"src/sprinx/data/cyto_cm/TRNAinf-{domain}-iso"
        records = []
        for aa in aas:
            stdout, stderr, rc = run(["cmfetch", cm_db, f"{domain}-{aa}"])
            if rc != 0:
                raise RuntimeError(f"cmfetch failed for {domain}-{aa}: {stderr.strip()}")
            with tempfile.NamedTemporaryFile("w", suffix=".cm") as fh:
                fh.write(stdout)
                fh.flush()
                seq = consensus_seq(fh.name)
                anticodon = derive_anticodon(seq, fh.name)
            records.append((aa, anticodon, seq))

        out_path = f"data/cyto/{domain}.fa"
        with open(out_path, "w") as fh:
            for aa, anticodon, seq in records:
                fh.write(f">{domain}_{aa}|{aa}|{anticodon}|synthetic_consensus\n{seq}\n")
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
