# sprinx

sprinx assigns Sprinzl coordinates to mitochondrial tRNA sequences. It is a
component of the qutrna2 pipeline, published as a standalone package for
transparency and reuse. It is not a general-purpose tRNA annotator and has
not been tested on cytosolic or bacterial tRNAs.

## What it does

Give sprinx a FASTA of mt-tRNA sequences (see "Header format" below for the
header conventions it understands) and it aligns each one to a covariance
model, works out whether an arm (D, T, or both) is missing or the
alignment just went wrong in that region, and reroutes truly armless
sequences to a matching armless model (Ozerova et al. 2024). It then assigns
a Sprinzl position to every nucleotide. The output is a per-nucleotide TSV:
Sprinzl position, structural region, which model was used, whether the
sequence was rerouted, and the arm-loss call. An optional standalone script
renders the result as an R2DT 2D structure plot.

## How it works

1. `cmalign` each sequence against every canonical CM tier
   (`--notrunc --nonbanded -g`). 
   Default: bacterial CM (endosymbiosis), then clade-specific CMs
   Keep the tier that anchors the anticodon and
   accounts for the most total base-paired columns across all stems (ties go
   to the earlier tier).
2. Locate the anticodon by position: the D-arm always precedes it, any
   variable arm and the T-arm always follow it. If it turns up in the first
   stem-loop instead of the second, the D-arm didn't occupy its own slot -
   D-arm missing.
3. Otherwise, check each stem: absent if it has fewer than 3 non-gap column
   pairs, or fewer than 2 of those are real WC/wobble pairs. This flags the
   D-arm, the T-arm, both, or neither.
4. For any arm flagged absent (other than a D-arm caught by step 2): does
   that span have enough sequence to physically close a hairpin (stem length
   + 3 nt)? If not, that's an arm loss. If so, fold just that span with
   RNAfold: a fold means cmalign mis-threaded the arm (patch it in place
   with the fold); no fold means arm loss after all.
5. Real arm loss: reroute to the matching armless CM (Ozerova et al.
   2024), isoacceptor ties broken by anticodon. Otherwise: assign Sprinzl
   labels on the alignment already in hand.

Exact thresholds and function names for each step are in the module
docstring at the top of [src/sprinx/label.py](src/sprinx/label.py).

## Installation

Install from a clone; `cmalign` (Infernal) needs to be on `PATH`
separately, see "Dependencies not on PyPI" below.

```bash
git clone <this repo>
cd sprinx
uv sync --extra viz              # or: pip install -e ".[viz]"
source .venv/bin/activate        # puts sprinx and this venv's python on PATH
```

The `--extra viz` / `[viz]` also installs `cairosvg`, needed only for
`scripts/visualize_ss.py`; leave it off if you only need the labeling TSV.
Activating the venv is what makes a plain `sprinx` and `python` (not your
system `python`) resolve correctly; `uv run <command>` is an equivalent
per-command alternative to activating.

## Quick start

Three commands, in order, from a FASTA of mt-tRNAs to a Sprinzl-labeled TSV,
QuTRNA2's `seq_to_sprinzl.tsv` format, and a PNG of the 2D structures. They
use the CMs and FASTA fixtures already checked into `data/`, so they run
as-is from a clone with no extra downloads (assumes the venv above is
activated, and `cmalign` and an R2DT Singularity image are on `PATH`; see
"Dependencies not on PyPI" below):

```bash
sprinx --fasta data/canonical.fa \
    --canonical-cm data/full_tRNAs_mitofinder_tRNAScanSE/TRNAinf-bact.cm \
                   data/full_tRNAs_mitofinder_tRNAScanSE \
    --armless-cm-dir data/truncated_cm/ \
    --out sprinzl_mapping.tsv

python scripts/convert_output_to_qutrna2-seq_to_sprinzl.py sprinzl_mapping.tsv

python scripts/visualize_ss.py --tsv sprinzl_mapping.tsv --out cloverleaves.png
```

## Dependencies not on PyPI

- `cmalign`, from [Infernal](http://eddylab.org/infernal/) >=1.1.4. Not a
  Python package; install via bioconda, a system package manager, or from
  source. (The `ViennaRNA` and `forgi` Python packages sprinx also needs are
  both on PyPI and install automatically with `pip install -e .`.)
- An [R2DT](https://r2dt.bio) Singularity image placed in `lib/r2dt`, only if you use
  `scripts/visualize_ss.py`.

## CM files (not bundled)

sprinx doesn't ship models for your organism; you supply them:

- **Canonical CMs**, via `--canonical-cm`: one or more directories of
  `{label}_{AA}.cm` files, or a single whole-family CM (e.g. Rfam RF00005, or
  a per-clade CM from MitoS2). Multiple sources are tried in order per
  sequence, first match wins; see "Why not just pick the best-scoring model?"
  below for why order, not score, decides.
- **Armless CMs**, via `--armless-cm-dir`: a directory of
  `armless_trn{AA}_wo_{d,t,d_and_t}.cm` files (Ozerova et al. 2024 naming).

Example, trying a bacterial whole-family model first (mitochondria's
bacterial ancestry makes it a good default guess), then a per-AA metazoan
directory:

```bash
sprinx --fasta my_mt_trnas.fa \
    --canonical-cm /path/to/cms/TRNAinf-bact.cm /path/to/cms/metazoan_per_aa/ \
    --armless-cm-dir /path/to/armless_cms/ \
    --out sprinzl_mapping.tsv
```

The repo's `data/` directory holds the CMs and FASTA fixtures the test suite
and the examples above use; see "Layout" below.

## Limitations

- The armless CMs (Ozerova et al. 2024) are mechanically truncated from
  canonical models, not retrained on armless sequences. They can mis-thread
  highly divergent armless mt-tRNAs.
- The order canonical CMs are tried in (e.g. bacterial whole-family, then
  metazoan per-AA, then armless rerouting) is a heuristic. It hasn't been
  tested outside metazoan mitochondrial sequences.
- The check that tells arm loss apart from a bad alignment works
  well for arm spans of roughly 13-20 nt; it hasn't been validated at the
  edges of that range.
- When exactly 3 stem-loops are found, the 2nd is always assumed to be the
  anticodon arm and the 3rd the T-arm. There's no way to instead read the
  3rd as a variable arm with the T-arm actually missing; that case isn't
  distinguished from an ordinary D-C-T cloverleaf with no variable arm.
- No support for cytosolic, bacterial, or archaeal tRNAs.

### Why not just pick the best-scoring model?

Because the score isn't comparable across models of different sizes. A
stripped-down armless model has fewer columns than a full canonical one, so
it scores canonical sequences better for reasons that have nothing to do with
biology, and normalizing by length doesn't fix this either, since an armless
model keeps the highest-information columns (the acceptor and anticodon
stems), which inflates its per-column score too. Picking a model by score or
E-value across models this different amounts to comparing numbers that were
never meant to be compared. sprinx instead tries one canonical model at a
time and only moves on when that model's alignment doesn't actually anchor
the anticodon. The full mechanism, including how a missing arm is told apart
from a misaligned one, is in the module docstring at the top of
[src/sprinx/label.py](src/sprinx/label.py); read it before touching the
arm-loss logic.

## Header format

Headers must use one of three forms:

- Pipe-delimited `id|aa|anticodon|taxon`, e.g. `seq1|Leu1|UAA|Mus_musculus`.
- An `anticodon=XXX` tag anywhere in the header, e.g. `seq1 anticodon=UAA`.
- GtRNAdb-style `tRNA-{AA}-{anticodon}` anywhere in the header, e.g.
  `mt-tRNA-Ala-TGC-1-1`.

The anticodon field is what drives model selection and arm-loss detection.
The aa field only picks which armless (or per-AA canonical) model family to
search; it has no other role. This matters for the GtRNAdb form specifically:
it never carries an isoacceptor digit, so Leu and Ser each cover two
anticodons under the same bare aa name. When that happens, sprinx tries each
matching model and keeps whichever one anchors the anticodon, the same
approach used to disambiguate filename-suffixed isoacceptor models
(Leu1/Leu2, Ser1/Ser2).

## Output

One row per nucleotide, written as TSV:

| column | meaning |
|---|---|
| `seq_id` | FASTA header |
| `seq_index` | 0-indexed position in the input sequence |
| `nucleotide` | base at that position |
| `sprinzl_position` | assigned label (`34`, `17a`, `60A`, ...) |
| `region` | structural region (`D_loop`, `T_stem_5`, `discriminator_CCA`, ...) |
| `cm_used` | which CM produced the final alignment |
| `rerouted` | whether the sequence got rerouted to an armless CM |
| `arm_loss_call` | structural diagnosis string; glossary below |
| `structure` | dot-bracket symbol at this position (sprinx's own final structure) |
| `cm_only_structure` | pre-patch structure at this position; blank unless this sequence needed an RNAfold patch |
| `rnafold_only_structure` | naive whole-sequence RNAfold structure at this position; blank unless this sequence needed an RNAfold patch |

The last three columns let
[scripts/visualize_ss.py](scripts/visualize_ss.py) render 2D diagrams
straight from the TSV, without re-running cmalign.

### Arm-loss call glossary

Every processed sequence gets exactly one of these:

- `CANONICAL_NO_ARM_LOSS`: every arm looks present.
- `T_OR_VAR_ARM_MISSING_slots=[n,..]`: one or more arm slots look empty
  (0-indexed, 5'->3'). A middle slot usually means an optional variable arm,
  not something sprinx reroutes for. The last slot is the T-arm: either
  truly missing, or patched via RNAfold if the alignment just misplaced
  it.
- `UPSTREAM_ARM_MISSING_offset=n`: the D-arm looks missing, caught by the
  anticodon landing further along the model than expected.
- `UPSTREAM_ARM_MISSING_slot=n`: the D-arm looks missing, but the anticodon
  didn't shift. Seen with CMs that model extra structure beyond the canonical
  D/C/T arms.
- `BOTH_ARMS_MISSING_slots=[n,..]`: both D-arm and T-arm look missing.
  Reroutes to `armless_trn{AA}_wo_d_and_t.cm`.
- `UNANCHORED_fallback_structurally_absent=[n,..]`: the anticodon couldn't be
  pinned down uniquely (an ambiguous AT-rich triplet), so this call is less
  reliable than the others.

A threading failure (alignment went wrong, arm isn't actually missing) is
logged as a separate line, not a call string: "CM diagnosed X-arm missing
(...) but the span folds as a real hairpin ... patching via RNAfold." The
patch is skipped, silently logged at DEBUG level, if it would conflict with
existing structure.

### Rendering 2D diagrams

Visualization is a separate standalone script, not part of the installable
package (R2DT needs a Singularity image, which is heavy and unnecessary for
anything just consuming sprinx's TSV output, e.g. QutRNA2):

```bash
python scripts/visualize_ss.py --tsv sprinzl_mapping.tsv --out cloverleaves.png
```

It draws one 2D diagram per sequence via R2DT, stitched into a single file:
`.svg`, `.png`, or `.pdf`, chosen by the extension on `--out` (R2DT itself
only emits SVG; PNG/PDF go through `cairosvg`). It plots sprinx's own final
structure per sequence, arm-loss calls and RNAfold patches included, rather
than a structure R2DT would work out on its own, which could disagree with
sprinx's diagnosis. For any sequence that got an RNAfold patch, two extra
files are also written, containing just those sequences: `_CMonly` (the
structure before the patch) and `_RNAfoldOnly` (the same sequence folded
naively as a whole, no CM at all), so the patch's effect is visible side by
side rather than assumed.

### Converting to QuTRNA2's format

```bash
python scripts/convert_output_to_qutrna2-seq_to_sprinzl.py sprinzl_mapping.tsv
# -> sprinzl_mapping.seq_to_sprinzl.tsv
```

Converts sprinx's output TSV into QuTRNA2's `seq_to_sprinzl.tsv` format: one
row per (Sprinzl label, tRNA id), giving that tRNA's 1-indexed sequence
position for the label, or `-` if the label doesn't occur in that sequence.
`id` is the FASTA header, unchanged.

## Layout

```
src/sprinx/
  label.py                   alignment, arm-loss classification, Sprinzl assignment
  cli.py                      argument parsing, per-record orchestration (console script `sprinx`)
scripts/
  visualize_ss.py             standalone R2DT 2D-diagram rendering, not part of the package
  convert_output_to_qutrna2-seq_to_sprinzl.py
                             converts sprinx's output TSV to QuTRNA2's seq_to_sprinzl.tsv format
recipe/
  meta.yaml                    conda recipe (bioconda-recipes conventions)
conftest.py                  pytest setup, loads .env / SPRINX_* vars for integration tests
env.example                  template for .env
data/                        example FASTA, canonical CM, armless CM library
  README.md                    data curation notes: evidence-tier definitions used in
                                curation_metadata.tsv
  curation_metadata.tsv        per-sequence literature evidence for each armless/doubly-armless
                                fixture entry (species, category, evidence tier, source DOI, notes)
  canonical.fa                 36 real cloverleaf mt-tRNAs (human, mouse, S. cerevisiae,
                                S. pombe), ground truth for "no arm-loss call should fire"
  D_armless.fa                 4 real D-armless mt-tRNAs, ground truth for missing_arm="d";
                                curated with literature evidence, see curation_metadata.tsv
  T_armless.fa                 4 real T-armless mt-tRNAs (Ascaris suum), ground truth for
                                missing_arm="t"; curated with literature evidence, see
                                curation_metadata.tsv
  both_armless.fa             3 real doubly-armless mt-tRNAs (R. culicivorax), ground truth
                                for missing_arm="d_and_t"; curated with literature evidence,
                                see curation_metadata.tsv
  spombe_mt.no_linker.fa     25 S. pombe mt-tRNAs (linker sequence trimmed), used as a
                                source of real sequences in integration tests
  TRNAinf-euk.cm                eukaryotic whole-family canonical CM (QutRNA2)
  mitofinder_models/           symlink to canonical CMs from MitoFinder, old INFERNAL-1 [1.0] format;
                                cmalign (Infernal 1.1.x) refuses these outright
  full_tRNAs_mitofinder_tRNAScanSE/ same CMs reformatted to current INFERNAL1/a via
                                `cmconvert -a`, one file in, one file out, same filename;
                                includes `TRNAinf-bact.cm`/`TRNAinf-euk.cm` (whole-family)
                                and per-AA `Metazoa_{AA}.cm` files; originals in
                                mitofinder_models/ are untouched; regenerate with:
                                `for f in data/mitofinder_models/*.cm; do
                                cmconvert -a "$f" > "data/full_tRNAs_mitofinder_tRNAScanSE/$(basename "$f")"; done`
  truncated_cm/                armless CM library, `armless_trn{AA}_wo_{d,t,d_and_t}.cm`
                                (Ozerova et al. 2024), used for --armless-cm-dir
  combined.cm*                 unused leftover from an earlier cmscan/combined-CM-database
                                exploration (see "Why not just pick the best-scoring model?"
                                above); not read by any current code path
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

Paths in `.env` must be absolute. Relative paths fail silently the moment
`cwd` differs from what you assumed.
