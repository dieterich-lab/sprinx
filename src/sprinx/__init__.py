"""
sprinx: Sprinzl-coordinate annotation for mitochondrial tRNAs.

assigns Sprinzl positions to mt-tRNA sequences by aligning each to a canonical
covariance model (cmalign), diagnosing which arm (D, T, or both) is genuinely
missing versus a CM threading failure via structural evidence rather than
alignment score, and rerouting genuinely armless sequences to the matching
armless CM (Ozerova et al. 2024) before numbering positions.

core logic: sprinx.label. CLI entry point: sprinx.cli:main (console script
`sprinx`). optional R2DT-rendered 2D diagrams of the output: the standalone
scripts/visualize_ss.py, not part of this package.
"""

__version__ = "0.1.0"
