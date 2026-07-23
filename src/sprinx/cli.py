"""
sprinx.cli: command-line entry point for sprinx.

argument parsing and per-record orchestration only; all labeling logic lives
in sprinx.label. optional R2DT-rendered visualization of the output TSV
is a separate standalone script, scripts/visualize_ss.py, not part of this
package (see its docstring for why: R2DT needs a Singularity image, which is
heavy and unnecessary for anything just consuming sprinx's TSV output).

usage: see README.md, or `sprinx --help`.
"""

import argparse
import multiprocessing
import os

import pandas as pd
from Bio import SeqIO
from loguru import logger

from sprinx.label import (
    _configure_logging,
    index_armless_cms,
    index_canonical_cms,
    process_one_record,
)


def main():
    parser = argparse.ArgumentParser(
        description="assign Sprinzl coordinates to mt-tRNA sequences via structure-based cm selection.")
    parser.add_argument("--fasta", required=True,
                        help="input FASTA; headers: 'id|aa|anticodon|taxon', 'anticodon=XXX' tag, "
                             "or GtRNAdb-style 'tRNA-{AA}-{anticodon}' name (e.g. mt-tRNA-Ala-TGC-1-1)")
    parser.add_argument("--canonical-cm", required=True, nargs="+", metavar="CM_OR_DIR",
                        help="one or more canonical CM sources, tried in order per sequence: "
                             "a path to a single CM (e.g. TRNAinf-bact.cm, applies to every aa), "
                             "or a directory of {label}_{AA}.cm files (e.g. Metazoan_P.cm) to "
                             "select per-sequence by aa. the first source whose alignment anchors "
                             "the anticodon unambiguously is used; earlier sources take priority "
                             "(e.g. a bacterial CM first, then a metazoan per-AA directory, since "
                             "a CM built for the wrong clade can fail to thread a divergent loop)")
    parser.add_argument("--armless-cm-dir", required=True,
                        help="directory (searched recursively) for "
                             "armless_trn{AA}_wo_{d,t,d_and_t}.cm files")
    parser.add_argument("--out", default="sprinzl_mapping.tsv",
                        help="output TSV path (default: sprinzl_mapping.tsv); includes a "
                             "'structure' column so scripts/visualize_ss.py can render it "
                             "without re-running cmalign")
    parser.add_argument("-p", "--processes", type=int, default=4,
                        help="worker processes (default: 4)")
    parser.add_argument("--debug", action="store_true",
                        help="log alignment, arm-loss diagnosis, and CM routing for every sequence")
    args = parser.parse_args()

    if args.debug:
        _configure_logging("DEBUG")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    records = [(str(r.id) + (" " + r.description.split(None, 1)[1]
                              if " " in r.description else ""),
                str(r.seq))
               for r in SeqIO.parse(args.fasta, "fasta")]
    armless_cm_index = index_armless_cms(args.armless_cm_dir)

    canonical_cm_tiers = []
    tier_descs = []
    for source in args.canonical_cm:
        if os.path.isdir(source):
            tier_index = index_canonical_cms(source)
            canonical_cm_tiers.append(tier_index)
            tier_descs.append(f"per-AA from {source} ({len(tier_index)} CMs)")
        else:
            canonical_cm_tiers.append(source)
            tier_descs.append(source)

    logger.info(f"{len(records)} sequences, canonical CM tiers (in priority order): "
                f"{tier_descs}, "
                f"{len(armless_cm_index)} armless CMs available for rerouting, "
                f"{args.processes} worker process(es)")

    tasks = [(header, seq, canonical_cm_tiers, armless_cm_index, args.debug)
             for header, seq in records]

    if args.processes > 1:
        with multiprocessing.Pool(args.processes) as pool:
            results = pool.map(process_one_record, tasks)
    else:
        # single-process path: --debug logging interleaves in real time without
        # multiprocess log-buffering surprises.
        results = [process_one_record(t) for t in tasks]

    all_rows = [row for r in results for row in r["rows"]]
    n_failed = sum(1 for r in results if not r["rows"])
    if n_failed:
        logger.warning(f"{n_failed}/{len(records)} sequences produced no output")

    pd.DataFrame(all_rows).to_csv(args.out, sep="\t", index=False)
    logger.info(f"table: {args.out}")


if __name__ == "__main__":
    main()
