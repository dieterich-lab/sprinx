"""
sprinx: Sprinzl-coordinate annotation for tRNAs.

for mt-tRNAs (--scheme mito): aligns each sequence to a canonical covariance
model (cmalign), diagnoses which arm (D, T, or both) is truly missing versus
a CM threading failure via structural evidence rather than alignment score,
and reroutes truly armless sequences to the matching armless CM (Ozerova et
al. 2024) before numbering positions.

structural parsing and Sprinzl-label assignment: sprinx.common (shared).
mito-specific CM tiering, arm-loss diagnosis, and armless-CM rerouting:
sprinx.mito. CLI entry point: sprinx.cli:main (console script `sprinx`).
optional R2DT-rendered 2D diagrams of the output: the standalone
scripts/visualize_ss.py, not part of this package.
"""

__version__ = "0.1.0"
