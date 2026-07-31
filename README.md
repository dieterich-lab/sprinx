# sprinx

Sprinzl numbering is the standard positional scheme for describing tRNA
structure across species. sprinx assigns Sprinzl coordinates to tRNA
sequences. It's a component of the QutRNA2 pipeline, usable standalone or
as a QutRNA2 dependency.

## What it does

Give sprinx a FASTA of tRNA sequences (see "Header format" below for the
header conventions it understands) and a `--scheme`:

- **`mito`**: aligns each sequence to a covariance model, works out whether
  an arm (D, T, or both) is missing or the alignment just went wrong in that
  region, and reroutes truly armless sequences to a matching armless model
  (Ozerova et al. 2024).
- **`euk` / `arch` / `bact`**: aligns each sequence to the matching
  per-isotype covariance model for that domain (from tRNAscan-SE). Cytosolic and
  nuclear tRNAs don't lose an arm the way mt-tRNAs do, so there's no
  arm-loss step on this path.

Either way, sprinx then assigns a Sprinzl position to every nucleotide. The
output is a per-nucleotide TSV: Sprinzl position, structural region, which
model was used, whether the sequence was rerouted, and the arm-loss call
(`mito` only). An optional standalone script renders the result as an R2DT
2D structure plot.

## Installation

Requires Python >=3.10. `cmalign` (Infernal) needs to be on `PATH`
separately, see "Dependencies not on PyPI" below.

From a clone:

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

From a pinned commit, without cloning (how qutrna2 consumes this today):

```bash
pip install git+https://github.com/dieterich-lab/sprinx.git@<commithash>
```

## Quick start

Three commands, in order, from a FASTA of mt-tRNAs to a Sprinzl-labeled TSV,
QuTRNA2's `seq_to_sprinzl.tsv` format, and a PNG of the 2D structures. The
default bundled CMs cover bacterial & metazoan mt-tRNAs, so no `--canonical-cm`/
`--armless-cm-dir` is needed here (see "CM files" below to override them).
This runs as-is from a clone with no extra downloads (assumes the venv above
is activated, `cmalign` is on `PATH`, and R2DT is reachable for the last
command; see "Dependencies not on PyPI" below):

```bash
sprinx --scheme mito --fasta data/mito/canonical.fa --out sprinzl_mapping.tsv

python scripts/convert_output_to_qutrna2-seq_to_sprinzl.py sprinzl_mapping.tsv

python scripts/visualize_ss.py --tsv sprinzl_mapping.tsv --out cloverleaves.png
```

## How it works

### `--scheme mito`

1. `cmalign` each sequence against every canonical CM tier
   (`--notrunc --nonbanded -g`; default tiers: a bacterial whole-family CM,
   then per-AA metazoan CMs). Keep the tier that anchors the anticodon and
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
docstring at the top of [src/sprinx/mito.py](src/sprinx/mito.py).

### `--scheme euk` / `arch` / `bact`

1. Resolve the header's aa field to the matching CM in that domain's
   combined per-isotype database (e.g. `euk-Ala`, `euk-iMet`, `euk-SeC`).
   `Ile2`/`iMet`/`fMet` are structurally distinct tRNAs that happen to
   share an amino acid identity, so this lookup matches the full field
   exactly.
2. `cmalign` directly against that one CM. Unlike the `mito` scheme,
   every cyto domain's database has one tRNAScan-SE CM per aa field,
   so there's no tier to fall back on if the alignment is poor.
3. Assign Sprinzl labels on the resulting alignment. Cytosolic/nuclear
   tRNAs don't lose arms the way mt-tRNAs do, so this path skips the
   arm-loss diagnosis entirely.

The CM naming convention and selection details are in the module docstring
at the top of [src/sprinx/cyto.py](src/sprinx/cyto.py).

### Stem re-seating (`--wc`)

A CM sometimes seats a helix one position off. Before labeling, sprinx checks
the neighbouring unpaired bases each internal stem could have paired with
instead, and slides the whole stem there if that gives more Watson-Crick or
wobble pairs. This is of concern for example when the stem is thermodynamically too short
to be stable, e.g. a 2-bp D-stem that may not need RNAfold to be patched into
a 3-bp stem using the adjacent bases.

`--wc 1` (default) checks the adjacent register, `2` also checks one position
further, `0` disables it. The anticodon stem never moves, since the numbering
is anchored to it. Changing this changes labels on D- and T-stem bases.

## Dependencies not on PyPI

- `cmalign`, from [Infernal](http://eddylab.org/infernal/) >=1.1.4. Not a
  Python package; install via bioconda, a system package manager, or from
  source. The `ViennaRNA` and `forgi` Python packages sprinx also needs are
  both on PyPI and install automatically with `pip install -e .`.
- [R2DT](https://r2dt.bio), only if you use `scripts/visualize_ss.py`. Either
  `r2dt.py` on `PATH`, or a container image run under docker, apptainer, or
  singularity. See [docs.r2dt.bio](https://docs.r2dt.bio) to obtain it.

## CM files

sprinx ships default CM databases, installed as package data so they're
available after a plain `pip install`:

- **`--scheme mito`** defaults to a bacterial whole-family CM first
  (mitochondria's bacterial ancestry makes it a good default guess), then a
  per-AA metazoan directory (both from MitoFinder/tRNAscan-SE). Covers
  metazoan mitochondrial tRNAs; a different clade needs its own CMs supplied
  via `--canonical-cm`/`--armless-cm-dir` (see below for the expected shape).
  See "Why not just pick the best-scoring model?" below for why order, not
  score, decides between multiple `--canonical-cm` sources.
  - For plant tRNAs (plastid or mitochondrial), add a eukaryotic CM as a
    middle tier, between the bacterial and metazoan defaults:
    ```bash
    sprinx --scheme mito --fasta plant_trnas.fa \
      --canonical-cm src/sprinx/data/mito_cm/canonical/TRNAinf-bact.cm \
                      data/mito/TRNAinf-euk.cm \
                      src/sprinx/data/mito_cm/canonical/mitofinder_models \
      --out sprinzl_mapping.tsv
    ```
- **`--scheme euk`/`arch`/`bact`** defaults to that domain's combined
  per-isotype CM database (tRNAscan-SE), one CM per amino acid. These aren't
  clade- or organism-specific, so no override is normally needed; use
  `--cyto-cm-db` only to point at a different database entirely.

To override the `mito` defaults for a different clade, supply your own
`--canonical-cm`/`--armless-cm-dir`; run `sprinx --help` for the exact
syntax each flag accepts.

The repo's `data/mito/` and `data/cyto/` directories hold the FASTA test
sequences the test suite and the examples above use (the CM databases
themselves live under `src/sprinx/data/`, since they ship as package data);
see "Layout" below.

## Limitations

- The `euk`/`arch`/`bact` schemes have no curated label set. CM selection and
  labeling are tested against GtRNAdb sequences (hundreds per domain:
  correct alignment, correct anticodon placement, structurally sound diagrams)
  and against each CM's consensus sequence. `euk` labels are tested
  against the nucleotides reported by Biela et al. 2023 review as conserved
  (`tests/data/conserved_positions.tsv`). That test reads the base
  under 11 positions. `arch` and `bact` have no such check.
  Cytosolic tRNAs are not expected to be truncated and are better conserved.
- The armless CMs (Ozerova et al. 2024) are mechanically truncated from
  canonical models, not retrained on armless sequences. They can mis-thread
  highly divergent armless mt-tRNAs.
- The order canonical CMs are tried in (e.g. bacterial whole-family, then
  metazoan per-AA, then armless rerouting) is a heuristic. It hasn't been
  tested outside metazoan mitochondrial sequences.
- The check that tells arm loss apart from a bad alignment folds the arm
  span with RNAfold. Below ~13 nt this hasn't been validated: competing
  folds may not be energetically negligible at that length, unlike the
  13-20 nt range this was checked against. The lower bound on span length
  itself (stem length + 3 nt) isn't a tunable threshold - the RNA backbone
  physically cannot close a shorter hairpin.
- When exactly 3 stem-loops are found, the 2nd is always assumed to be the
  anticodon arm and the 3rd the T-arm. There's no way to instead read the
  3rd as a variable arm with the T-arm missing; that case isn't
  distinguished from an ordinary D-C-T cloverleaf with no variable arm.
- A D-loop or T-loop shorter than expected drops positions in a fixed order,
  hardcoded in `common.py`. The order comes from the positions Biela et al. 2023
  reports as highly conserved. We checked it against Suzuki et al. 2020's 22
  curated human mt-tRNAs. It applies to every sequence, whatever the clade.
- An extra D-loop base takes 20a or 20b before 17a, because a eukaryotic D-loop
  grows at 20 rather than at 17. Placing it at 17a pushes the conserved G18 and
  G19 one slot later, which we measured on eukaryotic cytosolic tRNAs.
  - Only mito has some curated labels to check this against. Every cytosolic
    sequence we tested has a full 7-base T-loop, so T-loop shortening hasn't been
    tested for them.

### Why not just pick the best-scoring model?

Tried directly: `cmpress` every canonical and armless CM into one database,
`cmscan` each sequence, take the best E-value hit. On Ascaris and
Habronattus mt-tRNAs, wrong-isotype models outscored the correct one; for
Ascaris Asn (mtdbD00031155), the top two E-value hits are `armless_trnP_wo_t`
(Pro) and `H.seed25-1` (His); the true Asn model doesn't place in the top two
at all.

An armless CM has fewer columns than a canonical one, so it scores short
mt-tRNA sequences well for reasons unrelated to isotype match. E-value
corrects for database size only; it doesn't account for how much structure a
model attempts to score, so scores across differently-sized models aren't
comparable.

sprinx instead tries one canonical model at a time and moves on only when
the alignment fails to anchor the anticodon; see the module docstring at the
top of [src/sprinx/mito.py](src/sprinx/mito.py) for the full mechanism.

## Header format

Headers must use one of three forms:

- Pipe-delimited `id|aa|anticodon|taxon`, e.g. `seq1|Leu1|UAA|Mus_musculus`.
- An `anticodon=XXX` tag anywhere in the header, e.g. `seq1 anticodon=UAA`.
- GtRNAdb-style `tRNA-{AA}-{anticodon}` anywhere in the header, e.g.
  `mt-tRNA-Ala-TGC-1-1`.

The anticodon field is what drives model selection and arm-loss detection.
The aa field only picks which armless (or per-AA canonical) model family to
search; it has no other role. The GtRNAdb form never carries an isoacceptor
digit, so Leu and Ser each cover two anticodons under the same bare aa name.
When that happens, sprinx tries each matching model and keeps whichever one
anchors the anticodon, i.e. the same approach used to disambiguate
filename-suffixed isoacceptor models (Leu1/Leu2, Ser1/Ser2).

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

Every processed sequence gets exactly one call:

| call | meaning |
|---|---|
| `CANONICAL_NO_ARM_LOSS` | every arm looks present |
| `T_OR_VAR_ARM_MISSING_slots=[n,..]` | one or more arm slots look empty (0-indexed, 5'->3'); a middle slot is usually an optional variable arm; the last slot is the T-arm, missing or RNAfold-patched if the alignment misplaced it |
| `UPSTREAM_ARM_MISSING_offset=n` | D-arm missing, caught by the anticodon landing further along the model than expected |
| `UPSTREAM_ARM_MISSING_slot=n` | D-arm missing with no anticodon shift; seen with CMs modeling extra structure beyond D/C/T |
| `BOTH_ARMS_MISSING_slots=[n,..]` | D-arm and T-arm both missing; reroutes to `armless_trn{AA}_wo_d_and_t.cm` |
| `UNANCHORED_fallback_structurally_absent=[n,..]` | anticodon couldn't be pinned down uniquely (ambiguous AT-rich triplet); less reliable than the other calls |

A threading failure (alignment went wrong, arm isn't actually missing) logs
as a separate line: "CM diagnosed X-arm missing (...) but
the span folds as a real hairpin ... patching via RNAfold." A patch that
would conflict with existing structure is skipped and logged at DEBUG level.

### Rendering 2D diagrams

Visualization is a separate standalone script, not part of the installable
package. R2DT is a container, which is unnecessary for anything just
consuming sprinx's TSV output, e.g. QutRNA2. Example usage:

```bash
# docker on hand: nothing else to pass
python scripts/visualize_ss.py --tsv sprinzl_mapping.tsv --out cloverleaves.png

# a container image, e.g. on a cluster with apptainer or singularity
python scripts/visualize_ss.py --tsv sprinzl_mapping.tsv --out cloverleaves.png \
        --r2dt-image /path/to/r2dt.sif
```

`--r2dt-runtime` picks how `r2dt.py` is invoked: `auto` (the default) takes
`r2dt.py` from `PATH` if it is there, else the `--r2dt-image` container, else
docker with the `rnacentral/r2dt` image. Name a runtime explicitly with
`docker`, `apptainer`, `singularity`, or `native` to skip that search.

It draws one 2D diagram per sequence via R2DT, stitched into a single file:
`.svg`, `.png`, or `.pdf`, chosen by the extension on `--out` (R2DT itself
only emits SVG; PNG/PDF go through `cairosvg`). It plots sprinx's own final
structure per sequence, arm-loss calls and RNAfold patches included. This can
disagree with whatever structure R2DT would derive on its own from its
template library. For any sequence that got an RNAfold patch, two extra
files are also written, containing just those sequences: `_CMonly` (the
structure before the patch) and `_RNAfoldOnly` (the same sequence folded
naively as a whole, no CM at all), so the patch's effect is visible side by
side.

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
  common.py         shared structural parsing + Sprinzl-label assignment
  mito.py           canonical-CM tiering, arm-loss diagnosis, armless-CM rerouting
  cyto.py           combined-CM-database selection for --scheme euk/arch/bact
  cli.py            argument parsing, per-record orchestration (console script `sprinx`)
  data/             bundled package data
    mito_cm/          canonical + armless CMs for --scheme mito
    cyto_cm/          per-isotype CM databases for --scheme euk/arch/bact
scripts/
  visualize_ss.py                              standalone R2DT 2D-diagram rendering
  convert_output_to_qutrna2-seq_to_sprinzl.py   converts to QuTRNA2's seq_to_sprinzl.tsv format
  generate_synthetic_cyto_seqs.py               rebuilds data/cyto/{euk,arch,bact}.fa from the bundled CMs
  fetch_gtrnadb_seqs.py                         fetches real sequences into data/cyto/*_gtrnadb.fa
recipe/
  meta.yaml         conda recipe
conftest.py         pytest setup, loads .env / SPRINX_* vars for integration tests
env.example         template for .env
data/
  mito/             FASTA fixtures + curation notes for --scheme mito; see data/mito/README.md
  cyto/             test sequences for --scheme euk/arch/bact; see data/cyto/README.md
tests/
  test_sprinx_unit.py          unit tests, run anywhere
  test_sprinx_integration.py   runs real cmalign / RNAfold end to end
  data/
    test_data_bundle.txt       precomputed Stockholm alignments for the unit tests
    conserved_positions.tsv    nucleotides conserved across tRNAs, per Biela et al. 2023
output/             example run artifacts (TSVs, PNGs)
```

## Testing

`pytest` lives in the `test` extra.

```bash
uv run --extra test pytest tests/test_sprinx_unit.py
uv run --extra test pytest tests/test_sprinx_integration.py
                                             # requires cmalign and SPRINX_CANONICAL_CM /
                                             # SPRINX_ARMLESS_CM_DIR, set via .env
```

Paths in `.env` must be absolute. Relative paths fail silently the moment
`cwd` differs from what you assumed.

## License

GPL-3.0-or-later. See [LICENSE](LICENSE).
