#!/usr/bin/env python3
"""
visualize_ss.py: R2DT-rendered 2D diagrams for a sprinx TSV.

standalone script, not part of the installable sprinx package: R2DT is a
container, which is unnecessary for anything just consuming sprinx's TSV
output (e.g. QutRNA2). needs sprinx itself installed (for sprinx.common's
header parsing and subprocess helpers) plus its own extra dependency,
cairosvg (`pip install cairosvg`), and R2DT.

R2DT is not bundled. resolve_r2dt_runtime finds it as r2dt.py on PATH, or as a
container image run under apptainer, singularity, or docker. How to obtain it,
and the image name, mount point, and subcommand this script passes, are all
upstream's to define: see https://docs.r2dt.bio.

reads a sprinzl_mapping.tsv produced by `sprinx --out ...` (see sprinx.cli):
seq_id, seq_index, nucleotide, sprinzl_position, region, cm_used, rerouted,
arm_loss_call, structure, cm_only_structure, rnafold_only_structure. groups
rows back into per-record seq/ss/sprinzl, so no cmalign re-run is needed.

renders each record's final structure (the structure column) via R2DT's
template-free "stockholm" mode. this is sprinx's structural call, which
could disagree with whatever structure R2DT would derive on its own from its
template library. R2DT only accepts a real multi-sequence alignment as
input, but these records aren't aligned to each other at all, so
build_r2dt_stockholm fakes one: every record's sequence is concatenated
end-to-end into a single row, with one #=GC structureID region marking each
record's column span. see https://docs.r2dt.bio for the annotation
format.

for records where sprinx patched a CM threading failure via RNAfold (the
cm_only_structure / rnafold_only_structure columns are populated), also
renders the pre-patch CM-only structure and the naive whole-sequence RNAfold
fold alongside the primary plot, so the patch's effect is visible rather
than assumed.

usage: see README.md, or `python visualize_ss.py --help`.
"""

import argparse
import glob
import os
import shutil
import sys
import tempfile
import textwrap
import xml.etree.ElementTree as ET

import cairosvg
import pandas as pd
from loguru import logger

from sprinx.common import header_to_aa, header_to_taxon, run

# Upstream R2DT conventions, per https://docs.r2dt.bio: the published image
# name, and the mount point its docs bind the working directory to for both
# docker and singularity. Change here if upstream changes them.
R2DT_DEFAULT_DOCKER_IMAGE = "rnacentral/r2dt"
R2DT_CONTAINER_TEMP = "/rna/r2dt/temp"

R2DT_RUNTIMES = ("auto", "docker", "singularity", "apptainer", "native")

R2DT_SETUP_HINT = (
    "see https://docs.r2dt.bio to install R2DT, then either put r2dt.py on "
    "PATH or pass the container image to --r2dt-image"
)


class R2DTUnavailable(RuntimeError):
    """no usable R2DT runtime was found."""


def resolve_r2dt_runtime(runtime="auto", image=None):
    """pick how to invoke r2dt.py; returns (runtime_name, image_or_None).

    auto prefers a native r2dt.py on PATH, then a container image supplied via
    --r2dt-image, then docker with the public image. Raises R2DTUnavailable
    naming what was tried, rather than falling back to a path that may not
    exist."""
    if runtime == "native" or (runtime == "auto" and not image and shutil.which("r2dt.py")):
        if not shutil.which("r2dt.py"):
            raise R2DTUnavailable(f"r2dt.py is not on PATH.\n{R2DT_SETUP_HINT}")
        return "native", None

    if runtime in ("singularity", "apptainer"):
        if not image:
            raise R2DTUnavailable(f"--r2dt-runtime {runtime} needs --r2dt-image.\n{R2DT_SETUP_HINT}")
        if not shutil.which(runtime):
            raise R2DTUnavailable(f"{runtime} is not on PATH.\n{R2DT_SETUP_HINT}")
        return runtime, image

    if runtime == "docker":
        if not shutil.which("docker"):
            raise R2DTUnavailable(f"docker is not on PATH.\n{R2DT_SETUP_HINT}")
        return "docker", image or R2DT_DEFAULT_DOCKER_IMAGE

    # anything written as a path runs under apptainer/singularity; a bare name
    # like rnacentral/r2dt is a docker reference.
    if image and (os.path.isabs(image) or image.startswith((".", "~")) or os.path.isfile(image)):
        if not os.path.isfile(image):
            raise R2DTUnavailable(f"no R2DT image at {image}")
        for candidate in ("apptainer", "singularity"):
            if shutil.which(candidate):
                return candidate, image
        raise R2DTUnavailable(
            f"{image} needs apptainer or singularity, neither is on PATH.\n"
            f"{R2DT_SETUP_HINT}")
    if shutil.which("docker"):
        return "docker", image or R2DT_DEFAULT_DOCKER_IMAGE
    raise R2DTUnavailable(
        f"found no way to run R2DT: r2dt.py is not on PATH and docker is not "
        f"installed.\n{R2DT_SETUP_HINT}")


def r2dt_command(runtime, image, tmpdir, args):
    """full argv to run `r2dt.py <args>` under the chosen runtime. Container
    runtimes see tmpdir mounted at R2DT_CONTAINER_TEMP, so args must already
    use container-side paths; native mode substitutes the real tmpdir back in."""
    if runtime == "native":
        return ["r2dt.py"] + [a.replace(R2DT_CONTAINER_TEMP, tmpdir) for a in args]
    if runtime == "docker":
        # --user keeps the output owned by the caller, so the temp dir stays
        # removable; -w moves the working directory onto the bind mount, since
        # r2dt.py also writes paths relative to it.
        user = ["--user", f"{os.getuid()}:{os.getgid()}"] if hasattr(os, "getuid") else []
        return ["docker", "run", "--rm", *user, "-w", R2DT_CONTAINER_TEMP,
                "-v", f"{tmpdir}:{R2DT_CONTAINER_TEMP}", image, "r2dt.py"] + args
    return [runtime, "exec", "-B", f"{tmpdir}:{R2DT_CONTAINER_TEMP}", image,
            "r2dt.py"] + args


def _r2dt_id_line(segments):
    """build one #=GC structureID/regionID line: '|' marks each segment's first
    column, the segment's name fills the following columns, '.' pads the rest,
    the format R2DT's stockholm parser expects (see module docstring above)."""
    total = sum(length for _, length in segments)
    line = ["."] * total
    pos = 0
    for name, length in segments:
        line[pos] = "|"
        for i, ch in enumerate(name[:length - 1]):
            line[pos + 1 + i] = ch
        pos += length
    line[-1] = "|"
    return "".join(line)


def build_r2dt_stockholm(plotted):
    """records (each with 'header', 'seq', 'ss') -> (Stockholm text, region
    names) with one structureID region per record, in the same order as
    `plotted`. regionID is the aa field, so R2DT's --color-by region groups
    isoacceptors under one colour.

    names are a plain per-record index ('s0000', 's0001', ...), not derived
    from the header: a header-derived name must fit within that record's own
    sequence length (_r2dt_id_line's column span), so records sharing a long
    common prefix can truncate to the same name and make the later
    glob.glob(f"{name}_*.svg") in make_plot pick the wrong file. An index is
    always short enough to survive intact and unique by construction."""
    names = [f"s{i:04d}" for i in range(len(plotted))]

    concat_seq = "".join(r["seq"] for r in plotted)
    concat_ss = "".join(r["ss"] for r in plotted)
    structure_id = _r2dt_id_line([(name, len(r["seq"])) for name, r in zip(names, plotted)])
    region_id = _r2dt_id_line([(header_to_aa(r["header"]) or "na", len(r["seq"])) for r in plotted])

    sto_text = "\n".join([
        "# STOCKHOLM 1.0", "",
        f"seqs {concat_seq}", "",
        f"#=GC SS_cons     {concat_ss}",
        f"#=GC structureID {structure_id}",
        f"#=GC regionID    {region_id}",
        "//",
    ]) + "\n"
    return sto_text, names


def _cm_only_plot_path(path):
    """cloverleaves.svg -> cloverleaves_CMonly.svg, for the pre-RNAfold-patch
    comparison plot alongside the regular one."""
    root, ext = os.path.splitext(path)
    return f"{root}_CMonly{ext}"


def _rnafold_only_plot_path(path):
    """cloverleaves.svg -> cloverleaves_RNAfoldOnly.svg, for the naive
    whole-sequence-MFE comparison plot alongside the regular one."""
    root, ext = os.path.splitext(path)
    return f"{root}_RNAfoldOnly{ext}"


SVG_NS = "http://www.w3.org/2000/svg"
ET.register_namespace("", SVG_NS)


def _flip_panel_north(panel, width, height):
    """mirror a panel vertically in place (flip y only, x untouched).

    - why: R2R (the template-free layout engine R2DT uses here) always draws
      the acceptor stem at the bottom, with no orientation flag to change
      that (confirmed consistent across every sequence/shape checked, so one
      unconditional flip suffices).
    - why y-only, not a full 180-degree rotation: rotating would negate x
      too and swing every side arm from east to west. mirroring y alone
      moves the acceptor stem to the top while leaving east/west as R2R
      drew them.
    - every <text> glyph gets its own counter-mirror about its own y: two
      y-mirrors about different lines compose into a pure translation, so
      the glyph stays upright while still landing at its mirrored position -
      only the backbone/pairing geometry actually flips."""
    wrapper = ET.Element(f"{{{SVG_NS}}}g", {"transform": f"matrix(1 0 0 -1 0 {height})"})
    for child in list(panel):
        panel.remove(child)
        wrapper.append(child)
    panel.append(wrapper)
    for text in wrapper.iter(f"{{{SVG_NS}}}text"):
        ty = text.get("y", "0")
        text.set("transform", f"matrix(1 0 0 -1 0 {2 * float(ty)})")


NUCLEOTIDE_LETTERS = set("ACGUN")
ANTICODON_LABELS = {"34", "35", "36"}
CCA_TAIL_LABELS = {"74", "75", "76"}


def _is_block_boundary(seq_idx, region):
    """True if seq_idx is the first or last position of its region (its
    SPRINZL_REGION block), or has no region on one side (sequence end)."""
    prev_region = region.get(seq_idx - 1)
    next_region = region.get(seq_idx + 1)
    return region.get(seq_idx) != prev_region or region.get(seq_idx) != next_region


def _inject_sprinzl_labels(panel, sprinzl, region):
    """replace R2DT's own sequence-position numbering with sprinx's Sprinzl
    labels: shown at each SPRINZL_REGION block's start/end, every lettered
    insertion (17a, 20a, ...), and the anticodon (34/35/36) always. also
    colors the anticodon red and any present CCA-tail bases (74/75/76) gray.

    nucleotides are top-level <g><title>i (...)</title><text>BASE</text></g>
    in strict 5'->3' order (R2DT's own emission order); a running count of
    real base letters lines up with sprinzl's 0-indexed final_seq positions."""
    for g in list(panel):
        text = g.find(f"{{{SVG_NS}}}text")
        line = g.find(f"{{{SVG_NS}}}line")
        cls = (text.get("class") if text is not None else None) or \
              (line.get("class") if line is not None else None) or ""
        if "numbering-label" in cls or "numbering-line" in cls:
            panel.remove(g)

    seq_idx = 0
    for g in list(panel):
        text = g.find(f"{{{SVG_NS}}}text")
        if text is None or (text.text or "").strip() in ("5'", "3'", ""):
            continue
        if (text.text or "").strip().upper() not in NUCLEOTIDE_LETTERS:
            continue
        label = sprinzl.get(seq_idx, "")
        if label in ANTICODON_LABELS:
            text.set("style", "fill:red")
        elif label in CCA_TAIL_LABELS:
            text.set("style", "fill:gray")
        show = label and (not label.isdigit() or label in ANTICODON_LABELS
                          or _is_block_boundary(seq_idx, region))
        seq_idx += 1
        if not show:
            continue
        x, y = float(text.get("x")), float(text.get("y"))
        label_el = ET.SubElement(panel, f"{{{SVG_NS}}}text", {
            "x": str(x - 10), "y": str(y + 13),
            "class": "numbering-label",
        })
        label_el.text = label


CAPTION_FONT_SIZE = 11
CAPTION_LINE_HEIGHT = 14


def _wrap_caption(text, cell_w, max_lines=4):
    """text (header and summary, '\\n'-separated) -> wrapped lines that fit
    cell_w, each original line wrapped independently so the header and
    summary never run together into one blob. width is estimated from
    monospace glyph width since this is drawn as SVG <text>, not measured
    by a real layout engine."""
    chars_per_line = max(int(cell_w / (CAPTION_FONT_SIZE * 0.62)), 10)
    lines = [line for para in text.split("\n")
             for line in (textwrap.wrap(para, width=chars_per_line) or [""])]
    return lines[:max_lines]


def _grid_svg(panels, ncols, gap=20):
    """panels: list of (svg_root_element, width, height, caption, sprinzl,
    region).

    - arranges panels into a grid (ncols per row). R2DT's own --stitch
      only lays panels left-to-right in a single row, which gets unusably
      wide past a handful of sequences.
    - cell size is uniform (max panel width/height across all panels) so
      rows/columns stay aligned; each panel keeps its own native size,
      centered in its cell.
    - captions are wrapped to the cell width (see _wrap_caption) and
      centered above the panel: left-aligning to the cell would make a
      narrow panel in a wide cell read as belonging to its neighbour.
    - nested <svg> elements are SVG's own mechanism for embedding one
      diagram inside another at a given position/size; no rasterization
      needed to compose them."""
    cell_w = max(w for _, w, _, _, _, _ in panels) + gap
    caption_lines = [_wrap_caption(caption, cell_w) for _, _, _, caption, _, _ in panels]
    caption_height = max(len(lines) for lines in caption_lines) * CAPTION_LINE_HEIGHT + 6
    cell_h = max(h for _, _, h, _, _, _ in panels) + gap + caption_height
    nrows = -(-len(panels) // ncols)

    root = ET.Element(f"{{{SVG_NS}}}svg", {
        "width": str(ncols * cell_w), "height": str(nrows * cell_h),
    })
    for i, ((panel, w, h, _, sprinzl, region), lines) in enumerate(zip(panels, caption_lines)):
        row, col = divmod(i, ncols)
        x = col * cell_w + (cell_w - w) / 2
        y = row * cell_h + caption_height
        cx = x + w / 2
        for li, line in enumerate(lines):
            text = ET.SubElement(root, f"{{{SVG_NS}}}text", {
                "x": str(cx), "y": str(row * cell_h + 12 + li * CAPTION_LINE_HEIGHT),
                "font-family": "monospace", "font-size": str(CAPTION_FONT_SIZE),
                "text-anchor": "middle",
            })
            text.text = line
        _inject_sprinzl_labels(panel, sprinzl, region)
        _flip_panel_north(panel, w, h)
        panel.set("x", str(x))
        panel.set("y", str(y))
        root.append(panel)
    return ET.ElementTree(root)


def make_plot(records, out_path, runtime="auto", r2dt_image=None, ncols=6):
    """R2DT-rendered 2D diagram, one panel per record, arranged into our own
    grid (see _grid_svg: R2DT's own stitching doesn't wrap into rows).
    plotted in header order: species (taxon field) first, then tRNA (aa
    field), so isoacceptors of the same species group together and species
    cluster in the figure.

    Returns True when out_path was written. A caller that reports success
    must check it: R2DT failing, or emitting no SVG, leaves no file behind."""
    if not records:
        logger.warning("nothing to plot")
        return False
    plotted = sorted(records, key=lambda r: (header_to_taxon(r["header"]) or "",
                                              header_to_aa(r["header"]) or "", r["header"]))

    runtime, image = resolve_r2dt_runtime(runtime, r2dt_image)
    logger.debug(f"R2DT runtime: {runtime}" + (f", image {image}" if image else ""))

    sto_text, names = build_r2dt_stockholm(plotted)
    # a container can leave files the caller cannot delete; losing the temp dir
    # matters less than losing a finished plot.
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir:
        with open(os.path.join(tmpdir, "sprinx_plot.sto"), "w") as f:
            f.write(sto_text)
        cmd = r2dt_command(runtime, image, tmpdir,
                           ["stockholm", f"{R2DT_CONTAINER_TEMP}/sprinx_plot.sto",
                            f"{R2DT_CONTAINER_TEMP}/out", "--no-stitch"])
        _, stderr, rc = run(cmd)
        if rc != 0:
            logger.error(f"R2DT plotting failed: {stderr.strip()[:500]}")
            return False

        svg_dir = os.path.join(tmpdir, "out", "results", "svg")
        panels = []
        for name, r in zip(names, plotted):
            candidates = glob.glob(os.path.join(svg_dir, f"{name}_*.svg"))
            if not candidates:
                logger.warning(f"{r['header']}: R2DT produced no diagram for this sequence, skipping")
                continue
            root = ET.parse(candidates[0]).getroot()
            width, height = float(root.get("width")), float(root.get("height"))
            caption = f"{r['header'].split(None, 1)[0]}\n{r['summary']}"
            panels.append((root, width, height, caption, r["sprinzl"], r["region"]))

        if not panels:
            logger.error("R2DT plotting produced no SVG output")
            return False

        grid_path = os.path.join(tmpdir, "grid.svg")
        _grid_svg(panels, ncols).write(grid_path)

        ext = os.path.splitext(out_path)[1].lower()
        if ext == ".svg":
            shutil.copy(grid_path, out_path)
        else:
            _convert_svg(grid_path, out_path, ext)
        return True


_SVG_CONVERTERS = {".png": cairosvg.svg2png, ".pdf": cairosvg.svg2pdf}


# cairo's hard surface-size limit is ~32767px/side; a wide stitched plot (many
# sequences) can exceed that at scale=2.0, so the PNG scale is capped to keep
# the longer side under this, well clear of the real limit.
MAX_PNG_DIM = 16000


def _svg_intrinsic_size(svg_path):
    """(width, height) in px from an SVG's own width/height attributes."""
    root = ET.parse(svg_path).getroot()
    return float(root.get("width")), float(root.get("height"))


def _convert_svg(svg_path, out_path, ext):
    """R2DT only emits SVG; convert to whatever format --out asked for via
    cairosvg. .png and .pdf are supported; anything else falls back to .svg
    content copied as-is under the requested name, since silently writing
    nothing would be worse than an oddly-named SVG."""
    convert = _SVG_CONVERTERS.get(ext)
    if convert is None:
        logger.warning(f"--out: {ext} not supported for R2DT output (only .svg, "
                       f".png, .pdf); writing raw SVG to {out_path} instead")
        shutil.copy(svg_path, out_path)
        return
    kwargs = {}
    if ext == ".png":
        width, height = _svg_intrinsic_size(svg_path)
        kwargs = {"scale": min(2.0, MAX_PNG_DIM / max(width, height, 1))}
    try:
        convert(url=svg_path, write_to=out_path, **kwargs)
    except Exception as e:
        logger.warning(f"SVG->{ext} conversion failed: {e}")


def _records_from_tsv(tsv_path):
    """group a sprinx TSV's per-position rows back into per-record dicts:
    header, seq, ss, summary, sprinzl, region, cm_only_ss, rnafold_only_ss. a
    record with no rows in the TSV at all (sprinx skipped it upstream) is
    simply absent here, matching the original in-process filtering on empty
    rows."""
    df = pd.read_csv(tsv_path, sep="\t", keep_default_na=False)
    records = []
    for seq_id, group in df.groupby("seq_id", sort=False):
        group = group.sort_values("seq_index")
        cm_only_col = group["cm_only_structure"]
        rnafold_only_col = group["rnafold_only_structure"]
        cm_used = group["cm_used"].iloc[0]
        rerouted = bool(group["rerouted"].iloc[0])
        records.append({
            "header": seq_id,
            "seq": "".join(group["nucleotide"]),
            "ss": "".join(group["structure"]),
            "summary": f"CM:{cm_used}" + (" [rerouted]" if rerouted else ""),
            "sprinzl": {int(i): label for i, label in
                        zip(group["seq_index"], group["sprinzl_position"]) if label},
            "region": {int(i): r for i, r in
                       zip(group["seq_index"], group["region"])},
            "cm_only_ss": "".join(cm_only_col) if (cm_only_col != "").all() else None,
            "rnafold_only_ss": "".join(rnafold_only_col) if (rnafold_only_col != "").all() else None,
        })
    return records


def main():
    parser = argparse.ArgumentParser(
        description="render R2DT 2D diagrams from a sprinx TSV output.")
    parser.add_argument("--tsv", required=True,
                        help="sprinx TSV (sprinzl_mapping.tsv, from `sprinx --out ...`)")
    parser.add_argument("--out", required=True, metavar="SVG_OR_PNG_OR_PDF",
                        help="path for the R2DT-rendered 2D diagram, one panel per "
                             "sequence arranged in a grid; format is chosen by "
                             "extension (.svg, .png, .pdf)")
    parser.add_argument("--ncols", type=int, default=6, help="plot grid columns")
    parser.add_argument("--r2dt-image", metavar="IMAGE",
                        help="R2DT container image: a .sif path for apptainer or "
                             "singularity, or a docker image reference. defaults to "
                             f"{R2DT_DEFAULT_DOCKER_IMAGE} when docker is used")
    parser.add_argument("--r2dt-runtime", choices=R2DT_RUNTIMES, default="auto",
                        help="how to invoke r2dt.py. auto prefers r2dt.py on PATH, "
                             "then --r2dt-image, then docker (default: auto)")
    args = parser.parse_args()

    records = _records_from_tsv(args.tsv)
    try:
        ok = make_plot(records, args.out, runtime=args.r2dt_runtime,
                       r2dt_image=args.r2dt_image, ncols=args.ncols)
    except R2DTUnavailable as exc:
        logger.error(str(exc))
        sys.exit(1)
    if not ok:
        sys.exit(1)
    logger.info(f"plot: {args.out}")

    # sequences RNAfold-patched for a CM threading failure: also plot the
    # CM-only structure (pre-patch) side by side, so the patch's effect is
    # visible rather than assumed.
    cm_only_records = [{**r, "ss": r["cm_only_ss"]} for r in records if r.get("cm_only_ss")]
    if cm_only_records:
        cm_only_path = _cm_only_plot_path(args.out)
        if make_plot(cm_only_records, cm_only_path, runtime=args.r2dt_runtime,
                     r2dt_image=args.r2dt_image, ncols=args.ncols):
            logger.info(f"plot (CM-only, pre-RNAfold-patch): {cm_only_path}")

    # same sequences, but folded naively as a whole with RNAfold alone (no
    # CM at all): shows why the hybrid exists, since full-sequence MFE misses
    # tertiary contacts and modified bases a real mt-tRNA structure needs.
    rnafold_only_records = [{**r, "ss": r["rnafold_only_ss"]} for r in records if r.get("rnafold_only_ss")]
    if rnafold_only_records:
        rnafold_only_path = _rnafold_only_plot_path(args.out)
        if make_plot(rnafold_only_records, rnafold_only_path, runtime=args.r2dt_runtime,
                     r2dt_image=args.r2dt_image, ncols=args.ncols):
            logger.info(f"plot (RNAfold-only, whole-sequence naive fold): {rnafold_only_path}")


if __name__ == "__main__":
    main()
