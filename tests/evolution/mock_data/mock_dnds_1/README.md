# dN/dS Mock Data

This document explains the mock datasets, how they map to the script's inputs, how substitutions are reconstructed, and how the MEGA alignment corresponds to the script's view of ancestral → derived changes.

## Significant Sites (input to script)

The `significant_sites.tsv` file lists the sites considered significant by AlleleFlux. These rows are filtered by MAG, test type, and q-value in the script. Here is a preview:

| mag_id   | contig   |   position | gene_id   | test_type               |   min_p_value |   q_value |
|:---------|:---------|-----------:|:----------|:------------------------|--------------:|----------:|
| MAG_001  | contig1  |        161 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        108 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        185 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        102 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        110 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |         99 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        167 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        160 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        101 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        129 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |


## Profile Tables (pre/ancestral vs post/derived)


Both `pre_sample_MAG_001_profiled.tsv` and `post_sample_MAG_001_profiled.tsv` follow the same schema and provide base counts per (contig, position, gene_id). The key columns relevant to substitution reconstruction are:

- `contig`, `position`, `gene_id`: Site coordinates. `position` is 0-based on the contig.
- `ref_base`: Reference base in the forward-strand orientation of the contig.
- `A`, `T`, `G`, `C`: Observed counts used to pick the **major allele**.

**Major allele selection rule (ties):** choose the base with the highest count; on ties, prefer `ref_base`; if the ref_base is not in the tie, choose one among the tied alleles (script uses random choice; the mock data here doesn't rely on that randomness).


## Prodigal Genes

Prodigal headers encode 1-based start/end and strand. Sequences in the FASTA are **coding strand** (5'→3') regardless of genomic strand.

| gene_id   |   start |   end |   strand |   length_nt |
|:----------|--------:|------:|---------:|------------:|
| geneA     |     100 |   189 |        1 |          90 |
| geneB     |     600 |   671 |       -1 |          72 |



**Coordinate conversion used by the script** (`get_codon_from_site`):

- Forward gene: `pos_in_gene = position - (start - 1)`  
- Reverse gene: `pos_in_gene = (end - 1) - position`  

The codon containing a site starts at `codon_start_index = floor(pos_in_gene/3)*3`.

## Reconstructed Ancestral ORFs

Starting from the Prodigal coding sequences, the script substitutes the **pre** timepoint major allele at every profiled site. This yields the **ancestral ORF** per gene (still 5'→3' coding strand for both forward and reverse genes).

Concatenated ancestral ORF (geneA + geneB) matches the first sequence in MEGA exactly:


- Length: 162 nt
- Differences vs MEGA ancestral: 0


## MEGA Pairwise Alignment Consistency

`mega_pair_concat.fasta` contains two sequences (ancestral first, derived second). They correspond to the concatenation of `geneA` followed by `geneB` in coding orientation. We compared them to the script's reconstructed sequences:

- Differences (ancestral vs script): **0**
- Differences (derived vs script): **0**



## Site-Level Substitutions (forward-strand alleles)

For each significant site, we list the forward-strand major alleles at pre (ancestral) and post (derived) timepoints, then show the **ancestral and derived codons** after applying strand-aware logic (reverse-strand positions are complemented when incorporated into the coding sequence).

| gene_id   | contig   |   position |   strand |   pos_in_gene |   pos_in_codon | ancestral_major_fwd   | derived_major_fwd   | ancestral_codon   | derived_codon   | ancestral_aa   | derived_aa   | mutation_type   |
|:----------|:---------|-----------:|---------:|--------------:|---------------:|:----------------------|:--------------------|:------------------|:----------------|:---------------|:-------------|:----------------|
| geneA     | contig1  |        161 |        1 |            62 |              2 | T                     | C                   | TCT               | TCC             | S              | S            | S               |
| geneA     | contig1  |        108 |        1 |             9 |              0 | G                     | A                   | GCA               | ACA             | A              | T            | NS              |
| geneA     | contig1  |        185 |        1 |            86 |              2 | T                     | G                   | GCT               | GCG             | A              | A            | S               |
| geneA     | contig1  |        102 |        1 |             3 |              0 | G                     | C                   | GCC               | CCC             | A              | P            | NS              |
| geneA     | contig1  |        110 |        1 |            11 |              2 | A                     | T                   | GCA               | GCT             | A              | A            | S               |
| geneA     | contig1  |         99 |        1 |             0 |              0 | G                     | A                   | GCT               | ACT             | A              | T            | NS              |
| geneA     | contig1  |        167 |        1 |            68 |              2 | G                     | C                   | TCG               | TCC             | S              | S            | S               |
| geneA     | contig1  |        160 |        1 |            61 |              1 | C                     | A                   | TCT               | TAT             | S              | Y            | NS              |
| geneA     | contig1  |        101 |        1 |             2 |              2 | T                     | G                   | GCT               | GCG             | A              | A            | S               |
| geneA     | contig1  |        129 |        1 |            30 |              0 | C                     | G                   | CCG               | GCG             | P              | A            | NS              |
| geneA     | contig1  |        131 |        1 |            32 |              2 | G                     | T                   | CCG               | CCT             | P              | P            | S               |
| geneA     | contig1  |        157 |        1 |            58 |              1 | T                     | G                   | GTA               | GGA             | V              | G            | NS              |
| geneB     | contig1  |        635 |       -1 |            35 |              2 | A                     | G                   | GTT               | GTC             | V              | V            | S               |
| geneB     | contig1  |        604 |       -1 |            66 |              0 | G                     | A                   | CGA               | TGA             | R              | *            | NS              |
| geneB     | contig1  |        662 |       -1 |             8 |              2 | T                     | C                   | GGA               | GGG             | G              | G            | S               |
| geneB     | contig1  |        607 |       -1 |            63 |              0 | G                     | T                   | CGG               | AGG             | R              | R            | S               |
| geneB     | contig1  |        647 |       -1 |            23 |              2 | A                     | G                   | ACT               | ACC             | T              | T            | S               |
| geneB     | contig1  |        663 |       -1 |             7 |              1 | C                     | G                   | GGA               | GCA             | G              | A            | NS              |
| geneB     | contig1  |        659 |       -1 |            11 |              2 | A                     | T                   | CCT               | CCA             | P              | P            | S               |
| geneB     | contig1  |        645 |       -1 |            25 |              1 | G                     | A                   | ACC               | ATC             | T              | I            | NS              |


## Codon-Level Events (grouped)

The NG86 path-averaging in the script acts **per codon**. All significant sites that fall inside the same codon (same `gene_id` and `codon_start_index`) are combined, yielding `k` changes for that codon event. Below we show each codon event with the combined ancestral→derived codon and amino-acid change:

| gene_id   |   codon_start_index |   k | positions   | pos_in_codon   | ancestral_codon   | derived_codon   | ancestral_aa   | derived_aa   | mutation_type   |
|:----------|--------------------:|----:|:------------|:---------------|:------------------|:----------------|:---------------|:-------------|:----------------|
| geneA     |                   0 |   2 | 99,101      | 0,2            | GCT               | ACG             | A              | T            | NS              |
| geneA     |                   3 |   1 | 102         | 0              | GCC               | CCC             | A              | P            | NS              |
| geneA     |                   9 |   2 | 108,110     | 0,2            | GCA               | ACT             | A              | T            | NS              |
| geneA     |                  30 |   2 | 129,131     | 0,2            | CCG               | GCT             | P              | A            | NS              |
| geneA     |                  57 |   1 | 157         | 1              | GTA               | GGA             | V              | G            | NS              |
| geneA     |                  60 |   2 | 160,161     | 1,2            | TCT               | TAC             | S              | Y            | NS              |
| geneA     |                  66 |   1 | 167         | 2              | TCG               | TCC             | S              | S            | S               |
| geneA     |                  84 |   1 | 185         | 2              | GCT               | GCG             | A              | A            | S               |
| geneB     |                   6 |   2 | 663,662     | 1,2            | GGA               | GCG             | G              | A            | NS              |
| geneB     |                   9 |   1 | 659         | 2              | CCT               | CCA             | P              | P            | S               |
| geneB     |                  21 |   1 | 647         | 2              | ACT               | ACC             | T              | T            | S               |
| geneB     |                  24 |   1 | 645         | 1              | ACC               | ATC             | T              | I            | NS              |
| geneB     |                  33 |   1 | 635         | 2              | GTT               | GTC             | V              | V            | S               |
| geneB     |                  63 |   1 | 607         | 0              | CGG               | AGG             | R              | R            | S               |
| geneB     |                  66 |   1 | 604         | 0              | CGA               | TGA             | R              | *            | NS              |



### Notes on Synonymous / Non-synonymous calls
- The **mutation_type** above is based on translating the final ancestral vs derived codons using NCBI table 11 (bacterial/archaeal).
- For **k=1** codon events, this is identical to a site-level call.
- For **k>1**, the script computes NG86 **path-averaged** fractional S and NS counts by enumerating all orders of single-base changes and excluding paths with **intermediate stop codons**. The classification shown here (S or NS) simply reflects whether the *final* codon’s amino acid differs from the ancestral amino acid; the path-averaged S/NS **fractions** are used for dN/dS downstream.

## Edge Cases & Tie-Breaking


- **Ties in base counts:** broken in favor of the site’s `ref_base` from the profiles file.
- **Reverse-strand genes:** Profile alleles are in **forward-strand coordinates**, but the coding sequence for a reverse gene is the **reverse complement**; therefore, when composing codons, the script complements the forward-strand alleles before inserting into the codon.
- **Partial codons / out-of-bounds:** sites that would point outside the ORF or into partial codons are skipped with a warning.
- **Stops:** When path-averaging for **k>1**, any path that passes through an **intermediate stop** is excluded, but a final stop is allowed.

## Conclusion

The mock data and the script’s reconstruction are **consistent** with the MEGA pairwise alignment: the concatenated ancestral and derived sequences produced by the script **exactly match** the two sequences in `mega_pair_concat.fasta`. The site-level and codon-level substitution tables above show precisely which positions change, how they map into codons, and whether the resulting amino-acid change is synonymous or non-synonymous.
