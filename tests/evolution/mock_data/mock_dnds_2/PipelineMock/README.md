# Mock dN/dS Dataset mock_dnds_v4

This dataset contains **30 variable sites** across two coding genes (geneA forward, geneB reverse).
Sequences use 4-fold degenerate families to avoid introducing stop codons; events mix 3rd-position edits (often synonymous)
and 1st/2nd-position edits (often nonsynonymous). Some codons may have **k=2** edits to test NG86 path-averaging.

## MEGA files
- `MEGA_CDS/mega_pair_geneA.fasta`
- `MEGA_CDS/mega_pair_geneB.fasta`
- `MEGA_CDS/mega_pair_concat.fasta`

## Pipeline files (for your script)
- `prodigal_genes.fasta`
- `significant_sites.tsv` — 30 rows (0-based forward-strand coordinates on contig1)
- `profiles/pre_sample/pre_sample_MAG_001_profiled.tsv.gz`
- `profiles/post_sample/post_sample_MAG_001_profiled.tsv.gz`
