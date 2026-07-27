# GtRNAdb source data

Curated tRNA sequences fetched by `scripts/fetch_gtrnadb_seqs.py`, to
complement sprinx's synthetic-consensus cyto test data. Re-running the
script reproduces these files from the same URLs.

Sequences identical within a domain (common for multi-copy tRNA genes) are
collapsed to one record, keeping whichever header was encountered first.

## bact_gtrnadb.fa
99 unique sequences, 74 duplicate copies collapsed (see dedup_list.txt).
Sources:
- Bacillus_subtilis: https://gtrnadb.ucsc.edu/genomes/bacteria/Baci_subt_subtilis_168/baciSubt2-mature-tRNAs.fa
- Escherichia_coli: https://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/eschColi_K_12_MG1655-mature-tRNAs.fa

## arch_gtrnadb.fa
99 unique sequences, 15 duplicate copies collapsed (see dedup_list.txt).
Sources:
- Methanosarcina_barkeri: https://gtrnadb.ucsc.edu/genomes/archaea/Meth_bark_Fusaro/methBark1-mature-tRNAs.fa
- Haloferax_volcanii: https://gtrnadb.ucsc.edu/genomes/archaea/Halo_volc_DS2/haloVolc1-mature-tRNAs.fa

## euk_gtrnadb.fa
533 unique sequences, 767 duplicate copies collapsed (see dedup_list.txt).
Sources:
- Drosophila_melanogaster: https://gtrnadb.ucsc.edu/genomes/eukaryota/Dmela6/dm6-mature-tRNAs.fa
- Homo_sapiens: https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa
- Mus_musculus: https://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/mm39-mature-tRNAs.fa
- Schizosaccharomyces_pombe: https://gtrnadb.ucsc.edu/genomes/eukaryota/Schi_pomb_972h/schiPomb_972H-mature-tRNAs.fa
