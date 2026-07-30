"""
sprinx.cli: command-line entry point for sprinx.

argument parsing and per-record orchestration only; all labeling logic lives
in sprinx.common (shared) and sprinx.mito / sprinx.cyto (per-scheme). optional
R2DT-rendered visualization of the output TSV is a separate standalone script,
scripts/visualize_ss.py, not part of this package (see its docstring for why:
R2DT needs a Singularity image, which is heavy and unnecessary for anything
just consuming sprinx's TSV output).

--scheme selects which pipeline runs: mito uses tiered canonical-CM selection
with arm-loss diagnosis and armless-CM rerouting (sprinx.mito); euk/arch/bact
use combined-CM-database selection with no arm-loss step (sprinx.cyto), since
cytosolic/nuclear tRNAs don't lose arms the way mt-tRNAs do.

usage: see README.md, or `sprinx --help`.
"""

import argparse
import multiprocessing
import os

import pandas as pd
from Bio import SeqIO
from loguru import logger

from sprinx.common import _configure_logging, check_cm_source_formats

MITO_SCHEME = "mito"
CYTO_SCHEMES = ("euk", "arch", "bact")


def _load_records(fasta_path):
    return [(str(r.id) + (" " + r.description.split(None, 1)[1]
                          if " " in r.description else ""),
             str(r.seq))
            for r in SeqIO.parse(fasta_path, "fasta")]


def _run_pool(worker, tasks, processes):
    """dispatch tasks to worker, either via a process pool or in-process.
    single-process path exists so --debug logging interleaves in real time
    without multiprocess log-buffering surprises."""
    if processes > 1:
        with multiprocessing.Pool(processes) as pool:
            return pool.map(worker, tasks)
    return [worker(t) for t in tasks]


def _write_output(results, records, out_path):
    all_rows = [row for r in results for row in r["rows"]]
    n_failed = sum(1 for r in results if not r["rows"])
    if n_failed:
        logger.warning(f"{n_failed}/{len(records)} sequences produced no output")

    n_unanchored = sum(1 for r in results if r["rows"]
                       and r["rows"][0]["arm_loss_call"].startswith("UNANCHORED_fallback"))
    if n_unanchored:
        logger.warning(f"{n_unanchored}/{len(records)} sequences had an unanchored anticodon "
                       "(UNANCHORED_fallback) - less reliable than an anchored call, see README")

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    pd.DataFrame(all_rows).to_csv(out_path, sep="\t", index=False)
    logger.info(f"table: {out_path}")


def _run_mito(args, records):
    from sprinx.mito import (
        default_armless_cm_dir,
        default_canonical_cm_sources,
        index_armless_cms,
        index_canonical_cms,
        process_mito_record,
    )

    armless_cm_dir = args.armless_cm_dir or default_armless_cm_dir()
    canonical_cm_sources = args.canonical_cm or default_canonical_cm_sources()

    for source in [armless_cm_dir] + canonical_cm_sources:
        check_cm_source_formats(source)

    armless_cm_index = index_armless_cms(armless_cm_dir)

    canonical_cm_tiers = []
    tier_descs = []
    for source in canonical_cm_sources:
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

    tasks = [(header, seq, canonical_cm_tiers, armless_cm_index, args.debug, args.wc)
             for header, seq in records]
    return _run_pool(process_mito_record, tasks, args.processes)


def _run_cyto(args, records):
    from sprinx.cyto import default_cm_db_path, index_isotype_cms, process_cyto_record

    cm_db = args.cyto_cm_db or default_cm_db_path(args.scheme)
    check_cm_source_formats(cm_db)
    isotype_index = index_isotype_cms(cm_db)

    logger.info(f"{len(records)} sequences, isotype CM database: {cm_db} "
                f"({len(isotype_index)} CMs), {args.processes} worker process(es)")

    tasks = [(header, seq, cm_db, isotype_index, args.debug, args.wc)
             for header, seq in records]
    return _run_pool(process_cyto_record, tasks, args.processes)


def main():
    parser = argparse.ArgumentParser(
        description="assign Sprinzl coordinates to tRNA sequences via structure-based cm selection.")
    parser.add_argument("--scheme", required=True, choices=(MITO_SCHEME,) + CYTO_SCHEMES,
                        help="which pipeline to run: 'mito' (tiered canonical-CM "
                             "selection + arm-loss diagnosis + armless-CM rerouting), "
                             "or 'euk'/'arch'/'bact' (combined-CM-database selection, "
                             "no arm-loss step)")
    parser.add_argument("--fasta", required=True,
                        help="input FASTA; headers: 'id|aa|anticodon|taxon', 'anticodon=XXX' tag, "
                             "or GtRNAdb-style 'tRNA-{AA}-{anticodon}' name (e.g. mt-tRNA-Ala-TGC-1-1)")
    parser.add_argument("--canonical-cm", nargs="+", metavar="CM_OR_DIR",
                        help="(--scheme mito only) one or more canonical CM sources, tried in "
                             "order per sequence: a path to a single CM (e.g. TRNAinf-bact.cm, "
                             "applies to every aa), or a directory of {label}_{AA}.cm files "
                             "(e.g. Metazoan_P.cm) to select per-sequence by aa. the first "
                             "source whose alignment anchors the anticodon unambiguously is "
                             "used; earlier sources take priority (e.g. a bacterial CM first, "
                             "then a metazoan per-AA directory, since a CM built for the wrong "
                             "clade can fail to thread a divergent loop). default: sprinx's "
                             "bundled bacterial + metazoan per-AA CMs (covers metazoan "
                             "mitochondrial tRNAs; supply your own for a different clade)")
    parser.add_argument("--armless-cm-dir",
                        help="(--scheme mito only) directory (searched recursively) for "
                             "armless_trn{AA}_wo_{d,t,d_and_t}.cm files. default: sprinx's "
                             "bundled armless CM library (Ozerova et al. 2024)")
    parser.add_argument("--cyto-cm-db",
                        help="(--scheme euk/arch/bact only) path to that domain's pressed "
                             "combined CM database, one CM per amino acid. default: sprinx's "
                             "bundled tRNAscan-SE per-isotype database for the given --scheme")
    parser.add_argument("--out", default="sprinzl_mapping.tsv",
                        help="output TSV path (default: sprinzl_mapping.tsv); includes a "
                             "'structure' column so scripts/visualize_ss.py can render it "
                             "without re-running cmalign")
    parser.add_argument("-p", "--processes", type=int, default=4,
                        help="worker processes (default: 4)")
    parser.add_argument("--wc", type=int, choices=(0, 1, 2), default=1, metavar="N",
                        help="how far a stem may be re-seated when the same helix pairs "
                             "better elsewhere, correcting a CM that threaded it onto the "
                             "wrong bases. 0 disables, 1 (default) allows the adjacent "
                             "register, 2 also allows moving over a position. only applied "
                             "on a strict gain in Watson-Crick/wobble pairs, so a "
                             "well-threaded stem is left alone")
    parser.add_argument("--debug", action="store_true",
                        help="log alignment, arm-loss diagnosis, and CM routing for every sequence")
    args = parser.parse_args()

    if args.debug:
        _configure_logging("DEBUG")

    records = _load_records(args.fasta)

    if args.scheme == MITO_SCHEME:
        results = _run_mito(args, records)
    else:
        results = _run_cyto(args, records)

    _write_output(results, records, args.out)


if __name__ == "__main__":
    main()