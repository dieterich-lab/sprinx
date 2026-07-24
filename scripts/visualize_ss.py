#!/usr/bin/env python3
"""
visualize_ss.py: R2DT-rendered 2D diagrams for a sprinx TSV.

standalone script, not part of the installable sprinx package: R2DT needs a
Singularity image, which is unnecessary for anything just consuming
sprinx's TSV output (e.g. QutRNA2). needs sprinx itself installed (for
sprinx.common's header parsing and subprocess helpers) plus its own extra
dependency, cairosvg (`pip install cairosvg`), and a Singularity/R2DT image
(see README for setup).

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
import re
import shutil
import tempfile
import textwrap
import xml.etree.ElementTree as ET

import cairosvg
import pandas as pd
from loguru import logger

from sprinx.common import header_to_aa, header_to_taxon, run

R2DT_DEFAULT_IMAGE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "lib", "r2dt")


def _sanitize_r2dt_name(text, maxlen=40):
    """structureID/regionID names may not contain '|' or '.' (the annotation
    line's own delimiters); collapse anything else to '_' and truncate."""
    safe = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return safe[:maxlen] or "seq"


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
    isoacceptors under one colour. names are de-duplicated with a numeric
    suffix since sanitizing distinct headers can collide."""
    names, seen = [], {}
    for r in plotted:
        base = _sanitize_r2dt_name(r["header"])
        n = seen.get(base, 0)
        seen[base] = n + 1
        names.append(f"{base}_{n}" if n else base)

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
SPRINZL_LABEL_STEP = 5


def _inject_sprinzl_labels(panel, sprinzl, label_step=SPRINZL_LABEL_STEP):
    """replace R2DT's own plain sequence-position numbering (1, 2, 3, ...,
    shown every 10th residue by default, unrelated to Sprinzl coordinates)
    with sprinx's own Sprinzl labels.

    - shown every label_step-th integer position, but always for lettered
      insertions (17a, 20a, ...), since those don't follow a regular numeric
      cadence and would otherwise never appear.
    - each nucleotide is its own top-level <g><title>i (...)</title>
      <text>BASE</text></g>, emitted by R2DT in strict 5'->3' order. a
      running count of real base letters (skipping the synthetic 5'/3' end
      markers) lines up exactly with sprinzl's own 0-indexed final_seq
      positions."""
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
        seq_idx += 1
        show = label and (
            not label[:-1].isdigit()
            or int(re.match(r"\d+", label).group()) % label_step == 0
            or label == "1"
        )
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
    """panels: list of (svg_root_element, width, height, caption, sprinzl).

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
    cell_w = max(w for _, w, _, _, _ in panels) + gap
    caption_lines = [_wrap_caption(caption, cell_w) for _, _, _, caption, _ in panels]
    caption_height = max(len(lines) for lines in caption_lines) * CAPTION_LINE_HEIGHT + 6
    cell_h = max(h for _, _, h, _, _ in panels) + gap + caption_height
    nrows = -(-len(panels) // ncols)

    root = ET.Element(f"{{{SVG_NS}}}svg", {
        "width": str(ncols * cell_w), "height": str(nrows * cell_h),
    })
    for i, ((panel, w, h, _, sprinzl), lines) in enumerate(zip(panels, caption_lines)):
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
        _inject_sprinzl_labels(panel, sprinzl)
        _flip_panel_north(panel, w, h)
        panel.set("x", str(x))
        panel.set("y", str(y))
        root.append(panel)
    return ET.ElementTree(root)


def make_plot(records, out_path, r2dt_image=R2DT_DEFAULT_IMAGE, ncols=6):
    """R2DT-rendered 2D diagram, one panel per record, arranged into our own
    grid (see _grid_svg: R2DT's own stitching doesn't wrap into rows).
    plotted in header order: species (taxon field) first, then tRNA (aa
    field), so isoacceptors of the same species group together and species
    cluster in the figure. runs R2DT via its Singularity image (see README
    for setup); failures are logged and skipped, not raised, since this
    script is a sanity-check convenience, not sprinx's actual output."""
    if not records:
        logger.warning("nothing to plot")
        return
    plotted = sorted(records, key=lambda r: (header_to_taxon(r["header"]) or "",
                                              header_to_aa(r["header"]) or "", r["header"]))

    sto_text, names = build_r2dt_stockholm(plotted)
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(os.path.join(tmpdir, "sprinx_plot.sto"), "w") as f:
            f.write(sto_text)
        cmd = ["singularity", "exec", "-B", f"{tmpdir}:/rna/r2dt/temp", r2dt_image,
               "r2dt.py", "stockholm", "/rna/r2dt/temp/sprinx_plot.sto",
               "/rna/r2dt/temp/out", "--no-stitch"]
        _, stderr, rc = run(cmd)
        if rc != 0:
            logger.warning(f"R2DT plotting failed: {stderr.strip()[:500]}")
            return

        svg_dir = os.path.join(tmpdir, "out", "results", "svg")
        panels = []
        for name, r in zip(names, plotted):
            candidates = glob.glob(os.path.join(svg_dir, f"{name}_*.svg"))
            if not candidates:
                logger.warning(f"{r['header']}: R2DT produced no diagram for this sequence, skipping")
                continue
            root = ET.parse(candidates[0]).getroot()
            width, height = float(root.get("width")), float(root.get("height"))
            panels.append((root, width, height, f"{r['header']}\n{r['summary']}", r["sprinzl"]))

        if not panels:
            logger.warning("R2DT plotting produced no SVG output")
            return

        grid_path = os.path.join(tmpdir, "grid.svg")
        _grid_svg(panels, ncols).write(grid_path)

        ext = os.path.splitext(out_path)[1].lower()
        if ext == ".svg":
            shutil.copy(grid_path, out_path)
        else:
            _convert_svg(grid_path, out_path, ext)


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
    header, seq, ss, summary, sprinzl, cm_only_ss, rnafold_only_ss. a record
    with no rows in the TSV at all (sprinx skipped it upstream) is simply
    absent here, matching the original in-process filtering on empty rows."""
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
    parser.add_argument("--r2dt-image", default=R2DT_DEFAULT_IMAGE, metavar="PATH",
                        help=f"R2DT Singularity image (default: {R2DT_DEFAULT_IMAGE})")
    args = parser.parse_args()

    records = _records_from_tsv(args.tsv)
    make_plot(records, args.out, r2dt_image=args.r2dt_image, ncols=args.ncols)
    logger.info(f"plot: {args.out}")

    # sequences RNAfold-patched for a CM threading failure: also plot the
    # CM-only structure (pre-patch) side by side, so the patch's effect is
    # visible rather than assumed.
    cm_only_records = [{**r, "ss": r["cm_only_ss"]} for r in records if r.get("cm_only_ss")]
    if cm_only_records:
        cm_only_path = _cm_only_plot_path(args.out)
        make_plot(cm_only_records, cm_only_path, r2dt_image=args.r2dt_image, ncols=args.ncols)
        logger.info(f"plot (CM-only, pre-RNAfold-patch): {cm_only_path}")

    # same sequences, but folded naively as a whole with RNAfold alone (no
    # CM at all): shows why the hybrid exists, since full-sequence MFE misses
    # tertiary contacts and modified bases a real mt-tRNA structure needs.
    rnafold_only_records = [{**r, "ss": r["rnafold_only_ss"]} for r in records if r.get("rnafold_only_ss")]
    if rnafold_only_records:
        rnafold_only_path = _rnafold_only_plot_path(args.out)
        make_plot(rnafold_only_records, rnafold_only_path, r2dt_image=args.r2dt_image, ncols=args.ncols)
        logger.info(f"plot (RNAfold-only, whole-sequence naive fold): {rnafold_only_path}")


if __name__ == "__main__":
    main()
