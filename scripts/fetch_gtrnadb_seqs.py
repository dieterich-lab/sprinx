#!/usr/bin/env python3
"""
fetch_gtrnadb_seqs.py: pull real curated tRNA sequences from GtRNAdb, per
domain, to complement sprinx's synthetic-consensus cyto test data.

for each domain (euk/arch/bact), downloads one or more organisms' mature-tRNA
FASTA files, collapses identical sequences within a domain, and writes
data/cyto/{domain}_gtrnadb.fa. GtRNAdb's headers are left untouched -
header_to_aa/header_to_anticodon in common.py already parse that format.
Source URLs and collapsed duplicates go to data/cyto/README.md.

usage: python scripts/fetch_gtrnadb_seqs.py [--out-dir data/cyto]
"""
import argparse
import os
import re
import urllib.request

README_HEADER = """\
# GtRNAdb source data

Curated tRNA sequences fetched by `scripts/fetch_gtrnadb_seqs.py`, to
complement sprinx's synthetic-consensus cyto test data. Re-running the
script reproduces these files from the same URLs.

Sequences identical within a domain (common for multi-copy tRNA genes) are
collapsed to one record, keeping whichever header was encountered first.

Entries GtRNAdb tags 'Und' (undetermined isotype) or 'Sup' (suppressor
tRNA, reads a stop codon) are dropped: tRNAscan-SE's per-isotype CM
databases have no matching model for either.
"""

GTRNADB_SOURCES = {
    "bact": [
        ("Bacillus_subtilis",
         "https://gtrnadb.ucsc.edu/genomes/bacteria/Baci_subt_subtilis_168/baciSubt2-mature-tRNAs.fa"),
        ("Escherichia_coli",
         "https://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/eschColi_K_12_MG1655-mature-tRNAs.fa"),
    ],
    "arch": [
        ("Methanosarcina_barkeri",
         "https://gtrnadb.ucsc.edu/genomes/archaea/Meth_bark_Fusaro/methBark1-mature-tRNAs.fa"),
        ("Haloferax_volcanii",
         "https://gtrnadb.ucsc.edu/genomes/archaea/Halo_volc_DS2/haloVolc1-mature-tRNAs.fa"),
    ],
    "euk": [
        ("Drosophila_melanogaster",
         "https://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-mature-tRNAs.fa"),
        ("Homo_sapiens",
         "https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa"),
        ("Mus_musculus",
         "https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-mature-tRNAs.fa"),
        ("Schizosaccharomyces_pombe",
         "https://gtrnadb.ucsc.edu/genomes/eukaryota/Schi_pomb_972h/schiPomb_972H-mature-tRNAs.fa"),
    ],
}

# matches GtRNAdb's 'tRNA-{AA}-{anticodon}' substring, e.g. tRNA-Ala-GGC,
# tRNA-Ile2-CAT, tRNA-fMet-CAT, tRNA-Und-NNN (undetermined isotype call).
HEADER_RE = re.compile(r"tRNA-([A-Za-z]+\d*)-([A-Za-z]{3})")


def fetch_fasta(url):
    """download a FASTA file, return a list of (header, sequence) pairs.
    GtRNAdb rejects requests with no User-Agent header."""
    req = urllib.request.Request(url, headers={"User-Agent": "sprinx/0.1 (data fetch script)"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        text = resp.read().decode("utf-8")
    records, header, seq_lines = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_lines)))
            header, seq_lines = line[1:].strip(), []
        elif line.strip():
            seq_lines.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


# aa codes with no matching CM in tRNAscan-SE's per-isotype databases:
# 'Und' is GtRNAdb's undetermined-isotype call (anticodon 'NNN'); 'Sup'
# (suppressor tRNAs, real anticodon reading a stop codon) has no isotype
# bucket either, since amino-acid assignment is an anticodon->codon-table
# lookup and stop codons have no entry in that table.
NO_CM_AAS = {"Und", "Sup"}


def is_valid_record(header):
    """True if header names an aa with a matching CM and an ACGU anticodon.
    False for NO_CM_AAS. Raises on any other non-ACGU anticodon: an
    unexpected header format."""
    m = HEADER_RE.search(header)
    if not m:
        raise ValueError(f"header does not look like GtRNAdb format: {header!r}")
    aa, anticodon = m.group(1), m.group(2).upper().replace("T", "U")
    if aa in NO_CM_AAS:
        return False
    if re.fullmatch(r"[ACGU]{3}", anticodon):
        return True
    raise ValueError(f"unexpected non-ACGU anticodon in header: {header!r}")


def dedupe_by_sequence(records):
    """group records by sequence (case/T-U normalized). returns
    (deduped, duplicates_of): deduped is one (header, seq) per group, using
    the first-seen header; duplicates_of maps that header to the other
    headers collapsed into it (empty list if none)."""
    groups = {}
    for header, seq in records:
        key = seq.upper().replace("T", "U")
        groups.setdefault(key, {"seq": seq, "headers": []})["headers"].append(header)
    deduped = [(g["headers"][0], g["seq"]) for g in groups.values()]
    duplicates_of = {g["headers"][0]: g["headers"][1:] for g in groups.values()}
    return deduped, duplicates_of


def build_domain_fasta(domain, sources):
    records = []
    n_skipped = 0
    for _, url in sources:
        for header, seq in fetch_fasta(url):
            if is_valid_record(header):
                records.append((header, seq))
            else:
                n_skipped += 1
    if n_skipped:
        print(f"{domain}: skipped {n_skipped} entries with no isotype CM ({sorted(NO_CM_AAS)})")
    return dedupe_by_sequence(records)


def domain_readme_section(domain, sources, deduped, n_collapsed):
    lines = [
        f"\n## {domain}_gtrnadb.fa\n",
        (f"{len(deduped)} unique sequences, {n_collapsed} duplicate copies collapsed "
         f"(see dedup_list.txt).\n"),
        "Sources:\n",
    ]
    lines += [f"- {taxon}: {url}\n" for taxon, url in sources]
    return lines


def domain_dedup_lines(domain, duplicates_of):
    lines = [f"## {domain}_gtrnadb.fa\n"]
    dup_entries = [(h, d) for h, d in duplicates_of.items() if d]
    if not dup_entries:
        return lines + ["(no duplicates collapsed)\n"]
    for header, dups in dup_entries:
        lines.append(f"kept: {header}\n")
        lines += [f"  dup: {d}\n" for d in dups]
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", default="data/cyto")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    readme_sections = []
    dedup_sections = []
    for domain, sources in GTRNADB_SOURCES.items():
        deduped, duplicates_of = build_domain_fasta(domain, sources)
        out_path = os.path.join(args.out_dir, f"{domain}_gtrnadb.fa")
        with open(out_path, "w") as fh:
            fh.writelines(f">{header}\n{seq}\n" for header, seq in deduped)
        n_collapsed = sum(len(dups) for dups in duplicates_of.values())
        print(f"{domain}: {len(deduped)} unique sequences "
              f"({n_collapsed} duplicate copies collapsed) -> {out_path}")
        readme_sections += domain_readme_section(domain, sources, deduped, n_collapsed)
        dedup_sections += domain_dedup_lines(domain, duplicates_of)

    readme_path = os.path.join(args.out_dir, "README.md")
    with open(readme_path, "w") as fh:
        fh.write(README_HEADER)
        fh.writelines(readme_sections)
    print(f"wrote {readme_path}")

    dedup_path = os.path.join(args.out_dir, "dedup_list.txt")
    with open(dedup_path, "w") as fh:
        fh.writelines(dedup_sections)
    print(f"wrote {dedup_path}")


if __name__ == "__main__":
    main()