# sprinx

Assigns Sprinzl coordinates to mitochondrial tRNA sequences.

## Why this exists

mt-tRNAs fold into four different shapes: cloverleaf, D-armless, T-armless,
doubly-armless (Ozerova et al. 2024). Label positions relative to the
wrong shape and you get numbers that look plausible but are wrong. That's worse than
having no numbers at all.

The obvious approach is to align every sequence to every covariance model and keep
whichever alignment scores best. Don't do that. E-values are calibrated per model, and a
stripped-down armless CM has fewer columns than the canonical one, so it scores
canonical sequences better for reasons that have nothing to do with biology. Any
selection built on that comparison is built on sand.

Instead: align to the canonical CM once, read the actual base-pairing evidence, and make
a structural call about what's missing. Only reroute to a different model once you know
you need to.

## How it works

1. Align to the canonical CM with `cmalign --notrunc --nonbanded -g`. All three flags
   matter; see the module docstring before dropping one.
2. Anchor on the anticodon from the FASTA header, never inferred from position. A
   missing D-arm shifts the anticodon downstream in the model; a missing T-arm doesn't
   shift anything. That asymmetry tells you which arm is gone.
3. A stem slot with zero base-paired columns cannot form a stem, that's geometry, not a
   tunable threshold. But zero pairs has two causes, and they need different responses:
   the arm genuinely isn't there, or cmalign threaded a divergent sequence into insert
   columns instead of the model's stem columns. Counting nucleotides in the span against
   the physical minimum for a hairpin (stem length plus a 3nt loop) tells them apart.
4. Genuine arm loss reroutes to the matching `armless_trn{AA}_wo_{d,t,d_and_t}.cm`.
   Isoacceptors (Leu1/Leu2, Ser1/Ser2) get disambiguated by anticodon, not by filename.
5. Threading failures get patched, not rerouted: fold just the mis-threaded span
   (13-20nt) with RNAfold and splice the result in. Folding the whole molecule with
   RNAfold is a bad idea, mt-tRNA has tertiary contacts and modified bases that MFE
   doesn't account for.
6. Assign Sprinzl numbers across the resulting
   structure, including the D-armless replacement-loop case.

The full reasoning is in the module docstring at the top of
[sprinx.py](sprinx.py). Read it before touching the arm-loss logic.

## Output

One row per nucleotide, written as TSV:

| column | meaning |
|---|---|
| `seq_id` | FASTA header |
| `seq_index` | 0-indexed position in input sequence |
| `nucleotide` | base at that position |
| `sprinzl_position` | assigned label (`34`, `17a`, `60A`, ...) |
| `region` | structural region (`D_loop`, `T_stem_5`, `discriminator_CCA`, ...) |
| `cm_used` | which CM the final alignment came from |
| `rerouted` | whether it got rerouted to an armless CM |
| `arm_loss_call` | structural diagnosis string, glossary at the bottom of `sprinx.py` |

`--plot` renders one 2D diagram per sequence via [R2DT](https://r2dt.bio), stitched
into a single file -- `.svg`, `.png`, or `.pdf`, chosen by the extension on
`--plot`'s path (R2DT itself only emits SVG; PNG/PDF are converted from that via
`cairosvg`). It draws sprinx's own final structure for each sequence (arm-loss
calls and threading-failure patches included) using R2DT's template-free
`stockholm` mode -- not a structure R2DT would re-derive itself, which could
silently disagree with sprinx's own diagnosis. Since these sequences aren't a real
alignment, `build_r2dt_stockholm` fakes one: every sequence is concatenated
end-to-end into a single row, with one `#=GC structureID` region per sequence (see
[R2DT's Stockholm docs](https://docs.r2dt.bio/en/latest/stockholm-alignments.html)).
Useful for sanity-checking a run, not part of the actual output.

Requires an R2DT Singularity image (`--r2dt-image`, default `lib/r2dt` next to
`sprinx.py`) and `singularity`/`apptainer` on `PATH`.

For any sequence RNAfold-patched a CM threading failure, a second file is also
written with `_CMonly` inserted before the extension (e.g. `cloverleaves.svg` ->
`cloverleaves_CMonly.svg`), containing just those sequences rendered with their
pre-patch, CM-only structure -- so the patch's effect is visible side by side
rather than assumed.

## Requirements

- Python 3, with `pandas`, `RNA` (ViennaRNA), `forgi`, `biopython`, `loguru`,
  `scipy`, `cairosvg` (only needed for `--plot`'s PNG/PDF output).
- Infernal, with `cmalign` on `PATH`.
- For `--plot`: an R2DT Singularity image and `singularity`/`apptainer` on `PATH`.
- One or more canonical mt-tRNA CMs, e.g. `TRNAinf-euk.cm`. `--canonical-cm`
  takes multiple sources tried in priority order per sequence: the first
  source whose alignment anchors the anticodon unambiguously wins. Each
  source is either a single CM file (applies to every sequence, e.g. a
  whole-family CM like `TRNAinf-bact.cm`) or a directory of `{label}_{AA}.cm`
  files (e.g. `Metazoa_A.cm`), in which case the CM is chosen per-sequence by
  the header's aa field; `label` (clade or any prefix) is ignored. See
  `data/full_tRNAs_mitofinder_tRNAScanSE/`.

  This matters because a CM built for the wrong clade can lack the capacity
  to model a divergent loop (e.g. an unusually long variable loop before the
  T-stem); `cmalign` then threads the overflow into an adjacent arm's insert
  states and anticodon anchoring breaks entirely. Since mitochondria are of
  bacterial (endosymbiotic) origin, putting a bacterial whole-family CM ahead
  of a metazoan-only one can recover sequences the metazoan CM alone
  mis-threads -- selection is never by alignment score/E-value across tiers
  (see "Why this exists" above for why that's invalid), only by whether the
  anchor is clean.
- Armless CMs named `armless_trn{AA}_wo_{d,t,d_and_t}.cm`, see `data/truncated_cm/`.
  This naming convention comes from Ozerova et al. 2024; the indexer matches on it
  directly and skips anything that doesn't fit rather than guessing.

## Usage

```bash
python sprinx.py --fasta seqs.fa \
    --canonical-cm TRNAinf-euk.cm \
    --armless-cm-dir cm_models/ \
    --out results/sprinzl_mapping.tsv

# with plotting, parallel workers, and per-sequence decision logging
python sprinx.py --fasta seqs.fa \
    --canonical-cm TRNAinf-euk.cm \
    --armless-cm-dir cm_models/ \
    --plot results/cloverleaves.png \
    --processes 8 --debug

# --canonical-cm pointed at a directory: per-sequence canonical CM selection by aa
python sprinx.py --fasta data/canonical.fa \
    --canonical-cm data/full_tRNAs_mitofinder_tRNAScanSE \
    --armless-cm-dir data/truncated_cm/ \
    --out output/canonical_mitofinder.sprinzl.tsv \
    --plot output/canonical_mitofinder.png

# multiple --canonical-cm tiers, tried in order: bacterial whole-family CM first
# (mitochondria's endosymbiotic origin), falling back to the metazoan per-AA
# directory only for sequences the bacterial CM doesn't anchor cleanly
python sprinx.py --fasta data/canonical.fa \
    --canonical-cm data/full_tRNAs_mitofinder_tRNAScanSE/TRNAinf-bact.cm \
                   data/full_tRNAs_mitofinder_tRNAScanSE \
    --armless-cm-dir data/truncated_cm/ \
    --out output/canonical_mitofinder.sprinzl.tsv
```

Headers must be pipe-delimited as `id|aa|anticodon|taxon` (e.g.
`seq1|Leu1|UAA|Mus_musculus`), or carry an `anticodon=XXX` tag anywhere in the string.
The anticodon field drives CM selection and arm-loss anchoring. The aa field only
picks which armless CM family to search, it isn't load-bearing for anything structural.

## Layout

```
sprinx.py                   CLI, alignment, arm-loss classification, Sprinzl
                             assignment, plotting, all in one file
conftest.py                  pytest setup, loads .env / SPRINX_* vars for integration tests
env.example                  template for .env
data/                        example FASTA, canonical CM, armless CM library
  canonical.fa                 36 real cloverleaf mt-tRNAs (human, mouse, S. cerevisiae,
                                S. pombe), ground truth for "no arm-loss call should fire"
  D_armless.fa                 3 real D-armless mt-tRNAs, ground truth for missing_arm="d"
  T_armless.fa                 17 real T-armless mt-tRNAs, ground truth for missing_arm="t"
  both_armless_R_culicivorax_mt-tRNA-Ile.fa  1 real doubly-armless mt-tRNA (R. culicivorax),
                                ground truth for missing_arm="d_and_t"
  all.fa                       concatenation of the four files above (57 sequences),
                                for exercising every arm-loss shape in one run
  TRNAinf-euk.cm                eukaryotic whole-family canonical CM (QutRNA2)
  mitofinder_models/           symlink to canonical CMs from MitoFinder, old INFERNAL-1 [1.0] format;
                                cmalign (Infernal 1.1.x) refuses these outright
  full_tRNAs_mitofinder_tRNAScanSE/ same CMs reformatted to current INFERNAL1/a via
                                `cmconvert -a`, one file in, one file out, same filename;
                                includes `TRNAinf-bact.cm`/`TRNAinf-euk.cm` (whole-family)
                                and per-AA `Metazoa_{AA}.cm` files; originals in
                                mitofinder_models/ are untouched -- regenerate with:
                                `for f in data/mitofinder_models/*.cm; do
                                cmconvert -a "$f" > "data/full_tRNAs_mitofinder_tRNAScanSE/$(basename "$f")"; done`
  truncated_cm/                armless CM library, `armless_trn{AA}_wo_{d,t,d_and_t}.cm`
                                (Ozerova et al. 2024), used for --armless-cm-dir
tests/
  test_sprinx_unit.py          unit tests, run anywhere
  test_sprinx_integration.py   runs real cmalign / RNAfold end to end
  test_data_bundle.txt         precomputed Stockholm alignments for the unit tests
output/                      example run artifacts (TSVs, PNGs)
```

## Testing

```bash
pytest tests/test_sprinx_unit.py
pytest tests/test_sprinx_integration.py     # requires cmalign and SPRINX_CANONICAL_CM /
                                             # SPRINX_ARMLESS_CM_DIR, set via .env
```

Paths in `.env` must be absolute. Relative paths fail silently the moment `cwd` differs
from what you assumed.
