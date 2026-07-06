# sprinx

Assigns Sprinzl coordinates to mitochondrial tRNA sequences.

## Why this exists

mt-tRNAs fold into four different shapes: cloverleaf, D-armless, T-armless,
doubly-armless (Ozerova et al. 2024, PMC11571959). Label positions relative to the
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
6. Assign Sprinzl numbers (Sprinzl et al. 1998, PMC147216) across the resulting
   structure, including the D-armless replacement-loop case.

The full reasoning, with citations, is in the module docstring at the top of
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

`--plot` draws a grid of cloverleaves with Sprinzl labels, laid out with ViennaRNA's
NAVIEW. Useful for sanity-checking a run, not part of the actual output.

## Requirements

- Python 3, with `numpy`, `pandas`, `matplotlib`, `RNA` (ViennaRNA), `forgi`,
  `biopython`, `loguru`, `scipy`.
- Infernal, with `cmalign` on `PATH`.
- A canonical mt-tRNA CM, e.g. `TRNAinf-euk.cm`.
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
