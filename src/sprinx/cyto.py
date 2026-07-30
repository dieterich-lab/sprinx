"""
sprinx.cyto: combined-CM-database selection for cytosolic/nuclear tRNAs
(--scheme euk/arch/bact).

picks the matching per-isotype CM from a domain's combined database, then
hands the winning alignment straight to sprinx.common.sprinzl_map_from_alignment.

CM source: tRNAscan-SE's per-isotype combined covariance-model databases
(TRNAinf-{euk,arch,bact}-iso from
https://github.com/UCSC-LoweLab/tRNAscan-SE/tree/master/lib/models, checked
into data/cyto/isotype_cm/).

- one CM per amino acid field (e.g. euk-Ala, euk-Leu, euk-SeC, euk-iMet).
- Leu/Ser isoacceptors share a single CM per domain.
- Ile2 (AUA-decoding) and iMet/fMet (initiator) get their own model each:
  same amino acid, structurally distinct tRNA.

each database is one cmpress'd file holding every isotype's CM. cmalign
aligns to exactly one CM per call, so a candidate model is extracted into a
temp file via cmfetch first (see _cmfetch_one).

selection has no scoring step. the header's aa field picks exactly one CM
per domain (one model per amino acid). the anticodon anchor check still
runs, to warn on a mismatch; there's no fallback model to try instead.
"""

import os
import re
import tempfile

from loguru import logger

from sprinx.common import (
    _configure_logging,
    cmalign_one,
    finalize_structure,
    find_anticodon_stem_index,
    get_stem_loop_elements,
    header_to_aa,
    header_to_anticodon,
    package_data_path,
    run,
    slide_stems_to_improve_pairing,
    SPRINZL_REGION,
    sprinzl_map_from_alignment,
)

ISOTYPE_MODEL_RE = re.compile(r"^[a-z]+-([A-Za-z]+\d*)$")


def default_cm_db_path(domain):
    """bundled default --cyto-cm-db for a domain (euk/arch/bact): tRNAscan-SE's
    combined per-isotype CM database, one CM per amino acid."""
    return package_data_path("cyto_cm", f"TRNAinf-{domain}-iso")


def index_isotype_cms(cm_db_path):
    """scan a pressed tRNAscan-SE combined CM database for NAME records
    (one model per amino acid, e.g. 'euk-Ala') and return {aa_field: model_name}.

    aa_field is the raw suffix after the domain prefix: 'Ala', 'Ile2',
    'iMet', 'SeC'. matched directly against header_to_aa(header) - not
    passed through aa_field_to_cm_code's 3-letter conversion. these names
    already use the same short forms GtRNAdb-style headers use, and
    'iMet'/'fMet' aren't valid 3-letter amino acid codes at all."""
    stdout, stderr, rc = run(["grep", "^NAME", cm_db_path])
    if rc != 0:
        raise ValueError(f"couldn't read model names from {cm_db_path}: {stderr.strip()}")
    index = {}
    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[1]
        m = ISOTYPE_MODEL_RE.match(name)
        if not m:
            logger.debug(f"  not a recognized isotype model name, skipping: {name}")
            continue
        index[m.group(1)] = name
    logger.info(f"indexed {len(index)} isotype CMs in {cm_db_path}: {sorted(index)}")
    return index


def _cmfetch_one(cm_db_path, model_name):
    """extract one named model from a pressed multi-model CM database into a
    temp .cm file. cmalign needs a single-CM file; cmfetch/cmscan are the
    tools that address one model inside a larger pressed database."""
    stdout, stderr, rc = run(["cmfetch", cm_db_path, model_name])
    if rc != 0:
        logger.warning(f"cmfetch failed (rc={rc}) for {model_name} in {cm_db_path}:\n{stderr.strip()}")
        return None
    with tempfile.NamedTemporaryFile("w", suffix=".cm", delete=False) as fh:
        fh.write(stdout)
        return fh.name


def select_cyto_cm_and_align(header, seq, cm_db_path, isotype_index):
    """Top-level CM selection for one cytosolic/nuclear tRNA sequence.

    1. Resolve the header's aa field to a model name via isotype_index.
       Exact match only: Ile/Ile2 and Met/iMet/fMet are distinct tRNAs
       sharing an amino acid identity, so no isoacceptor-digit stripping.
    2. Align against that one model. A domain's combined database has
       exactly one CM per aa field: nothing left to choose between.
    3. aa field not in the index, or cmfetch/cmalign fails: log a warning,
       return no alignment.

    Returns: dict with final_alignment, cm_used (the model name, e.g.
    'euk-Ala' - there's no per-model file on disk to point to)."""
    aa_field = header_to_aa(header)
    model_name = isotype_index.get(aa_field)
    if model_name is None:
        logger.warning(f"{header}: aa field {aa_field!r} not found in isotype index "
                       f"{sorted(isotype_index)}; skipping")
        return {"final_alignment": None, "cm_used": None}

    cm_path = _cmfetch_one(cm_db_path, model_name)
    if cm_path is None:
        return {"final_alignment": None, "cm_used": model_name}
    try:
        aln = cmalign_one(header, seq, cm_path)
    finally:
        os.unlink(cm_path)

    if aln is None:
        return {"final_alignment": None, "cm_used": model_name}

    anticodon = header_to_anticodon(header)
    elements = get_stem_loop_elements(aln["ss_cons"])
    idx, method = find_anticodon_stem_index(aln["aligned_seq"], elements, anticodon)
    if idx is None:
        logger.warning(f"{header}: anticodon did not anchor against {model_name} "
                       f"({method}); using the alignment anyway, no fallback model available\n"
                       f"  aligned_seq={aln['aligned_seq']}\n"
                       f"  ss_cons={aln['ss_cons']}")

    return {"final_alignment": aln, "cm_used": model_name}


def process_cyto_record(args):
    """worker for one (header, seq) FASTA record, cytosolic/nuclear path.

    takes a single tuple for Pool.map compatibility. the winning alignment
    goes straight to sprinzl_map_from_alignment, with no arm-loss step of any kind."""
    header, seq, cm_db_path, isotype_index, debug, wc = args
    seq = seq.upper().replace("T", "U")

    if debug:
        _configure_logging("DEBUG")

    routing = select_cyto_cm_and_align(header, seq, cm_db_path, isotype_index)
    alignment = routing["final_alignment"]

    if alignment is None:
        logger.warning(f"{header}: cmalign failed, skipped")
        return {"header": header, "rows": [], "summary": "CMALIGN_FAILED"}

    final_seq, final_ss = finalize_structure(alignment)

    if len(final_seq) != len(seq):
        logger.warning(f"{header}: ungapped length {len(final_seq)} != input {len(seq)}, skipped")
        return {"header": header, "rows": [], "summary": "LENGTH_MISMATCH"}

    anticodon = header_to_anticodon(header)
    if anticodon is None:
        logger.warning(f"{header}: no anticodon in header; C-stem location unreliable")

    sprinzl = sprinzl_map_from_alignment(alignment, anticodon, missing_arm=None, wc=wc)

    unlabeled = [i for i in range(len(final_seq)) if i not in sprinzl]
    if unlabeled:
        logger.warning(f"{header}: {len(unlabeled)} position(s) left without a Sprinzl "
                       f"number at seq index {unlabeled}; output rows for them are blank")

    cm_name = routing["cm_used"] or "NONE"
    summary = f"CM:{cm_name}"

    logger.debug(f"{header}")
    logger.debug(f"  seq ({len(final_seq)}nt): {final_seq}")
    logger.debug(f"  ss  ({len(final_ss)}nt):  {final_ss}")
    logger.debug(f"  raw stockholm:\n{alignment['raw_sto']}")

    rows = []
    for i, base in enumerate(final_seq):
        label = sprinzl.get(i, "")
        region_key = re.match(r"e?\d+", label).group() if label else ""
        rows.append({
            "seq_id": header, "seq_index": i, "nucleotide": base,
            "sprinzl_position": label, "region": SPRINZL_REGION.get(region_key, ""),
            "cm_used": cm_name, "rerouted": False, "arm_loss_call": "",
            "structure": final_ss[i],
            "cm_only_structure": "", "rnafold_only_structure": "",
        })

    logger.info(f"{header}: {summary}")
    return {"header": header, "rows": rows, "summary": summary,
            "seq": final_seq, "ss": final_ss, "sprinzl": sprinzl}
