"""__init__.py - Sprinzl-coordinate annotation for tRNAs

input:   tRNA sequences in FASTA, cytosolic or mitochondrial
output:  a per-position Sprinzl table, written by the CLI
usage:   import sprinx    (console script: sprinx)
env:     cmalign, cmstat and cmfetch on PATH
notes:   common.py parses structure and assigns labels; mito.py handles CM
         tiering, arm-loss diagnosis and armless rerouting; cyto.py picks a
         per-isotype CM; cli.py is the entry point. scripts/visualize_ss.py
         renders 2D diagrams and is outside the package
"""

__version__ = "0.2.0"
