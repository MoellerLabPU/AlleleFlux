# dN/dS Mock Data

## Significant Sites (input to the script)

| mag_id   | contig   |   position | gene_id   | test_type               |   min_p_value |   q_value |
|:---------|:---------|-----------:|:----------|:------------------------|--------------:|----------:|
| MAG_001  | contig1  |        152 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        222 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        110 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        159 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        113 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        207 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        161 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        174 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        224 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |
| MAG_001  | contig1  |        151 | geneA     | two_sample_paired_tTest |         1e-06 |     1e-06 |


## Profile Tables (pre/ancestral vs post/derived)

Columns used: `contig`, `position` (0-based), `gene_id`, `ref_base`, and base counts `A/T/G/C`.

**Pre (ancestral) preview:**

| contig   |   position | gene_id   | ref_base   |   A |   T |   G |   C |
|:---------|-----------:|:----------|:-----------|----:|----:|----:|----:|
| contig1  |        152 | geneA     | C          |   0 |   0 |   0 |  30 |
| contig1  |        222 | geneA     | A          |  30 |   0 |   0 |   0 |
| contig1  |        110 | geneA     | A          |  30 |   0 |   0 |   0 |
| contig1  |        159 | geneA     | T          |   0 |  30 |   0 |   0 |
| contig1  |        113 | geneA     | T          |   0 |  30 |   0 |   0 |
| contig1  |        207 | geneA     | C          |   0 |   0 |   0 |  30 |
| contig1  |        161 | geneA     | T          |   0 |  30 |   0 |   0 |
| contig1  |        174 | geneA     | C          |   0 |   0 |   0 |  30 |
| contig1  |        224 | geneA     | C          |   0 |   0 |   0 |  30 |
| contig1  |        151 | geneA     | T          |   0 |  30 |   0 |   0 |

**Post (derived) preview:**

| contig   |   position | gene_id   | ref_base   |   A |   T |   G |   C |
|:---------|-----------:|:----------|:-----------|----:|----:|----:|----:|
| contig1  |        152 | geneA     | C          |   0 |  25 |   0 |   5 |
| contig1  |        222 | geneA     | A          |   5 |  25 |   0 |   0 |
| contig1  |        110 | geneA     | A          |   5 |   0 |  25 |   0 |
| contig1  |        159 | geneA     | T          |   0 |   5 |  25 |   0 |
| contig1  |        113 | geneA     | T          |   0 |   5 |   0 |  25 |
| contig1  |        207 | geneA     | C          |  25 |   0 |   0 |   5 |
| contig1  |        161 | geneA     | T          |   0 |   5 |   0 |  25 |
| contig1  |        174 | geneA     | C          |   0 |   0 |  25 |   5 |
| contig1  |        224 | geneA     | C          |  25 |   0 |   0 |   5 |
| contig1  |        151 | geneA     | T          |   0 |   5 |   0 |  25 |


## Prodigal Genes

Prodigal FASTA headers encode 1-based start, end, and strand. Sequences are coding-strand (5'→3').

| gene_id   |   start |   end |   strand |   length_nt |
|:----------|--------:|------:|---------:|------------:|
| geneA     |     100 |   225 |        1 |         126 |
| geneB     |     900 |  1007 |       -1 |         108 |


## MEGA Pairwise Alignment Consistency

- Ancestral concat length: **234**; MEGA ancestral length: **234**; differences: **0**

- Derived concat length: **234**; MEGA derived length: **234**; differences: **0**



## Site-Level Substitutions (forward-strand alleles)

| gene_id   | contig   |   position |   strand |   pos_in_gene |   pos_in_codon | ancestral_major_fwd   | derived_major_fwd   | ancestral_codon   | derived_codon   | ancestral_aa   | derived_aa   | mutation_type   |
|:----------|:---------|-----------:|---------:|--------------:|---------------:|:----------------------|:--------------------|:------------------|:----------------|:---------------|:-------------|:----------------|
| geneA     | contig1  |        110 |        1 |            11 |              2 | A                     | G                   | GCA               | GCG             | A              | A            | S               |
| geneA     | contig1  |        113 |        1 |            14 |              2 | T                     | C                   | GGT               | GGC             | G              | G            | S               |
| geneA     | contig1  |        116 |        1 |            17 |              2 | C                     | T                   | GGC               | GGT             | G              | G            | S               |
| geneA     | contig1  |        120 |        1 |            21 |              0 | G                     | A                   | GGA               | AGA             | G              | R            | NS              |
| geneA     | contig1  |        124 |        1 |            25 |              1 | C                     | T                   | CCT               | CTT             | P              | L            | NS              |
| geneA     | contig1  |        126 |        1 |            27 |              0 | C                     | G                   | CCC               | GCC             | P              | A            | NS              |
| geneA     | contig1  |        151 |        1 |            52 |              1 | T                     | C                   | GTC               | GCC             | V              | A            | NS              |
| geneA     | contig1  |        152 |        1 |            53 |              2 | C                     | T                   | GTC               | GTT             | V              | V            | S               |
| geneA     | contig1  |        159 |        1 |            60 |              0 | T                     | G                   | TCT               | GCT             | S              | A            | NS              |
| geneA     | contig1  |        161 |        1 |            62 |              2 | T                     | C                   | TCT               | TCC             | S              | S            | S               |
| geneA     | contig1  |        173 |        1 |            74 |              2 | T                     | G                   | CGT               | CGG             | R              | R            | S               |
| geneA     | contig1  |        174 |        1 |            75 |              0 | C                     | G                   | CGC               | GGC             | R              | G            | NS              |
| geneA     | contig1  |        182 |        1 |            83 |              2 | A                     | C                   | CGA               | CGC             | R              | R            | S               |
| geneA     | contig1  |        187 |        1 |            88 |              1 | C                     | A                   | GCC               | GAC             | A              | D            | NS              |
| geneA     | contig1  |        206 |        1 |           107 |              2 | A                     | T                   | GGA               | GGT             | G              | G            | S               |
| geneA     | contig1  |        207 |        1 |           108 |              0 | C                     | A                   | CCT               | ACT             | P              | T            | NS              |
| geneA     | contig1  |        222 |        1 |           123 |              0 | A                     | T                   | ACC               | TCC             | T              | S            | NS              |
| geneA     | contig1  |        224 |        1 |           125 |              2 | C                     | A                   | ACC               | ACA             | T              | T            | S               |
| geneB     | contig1  |        994 |       -1 |            12 |              0 | G                     | C                   | CCC               | GCC             | P              | A            | NS              |
| geneB     | contig1  |        989 |       -1 |            17 |              2 | C                     | G                   | CCG               | CCC             | P              | P            | S               |
| geneB     | contig1  |        977 |       -1 |            29 |              2 | C                     | T                   | ACG               | ACA             | T              | T            | S               |
| geneB     | contig1  |        974 |       -1 |            32 |              2 | T                     | C                   | ACA               | ACG             | T              | T            | S               |
| geneB     | contig1  |        949 |       -1 |            57 |              0 | G                     | A                   | CGT               | TGT             | R              | C            | NS              |
| geneB     | contig1  |        947 |       -1 |            59 |              2 | A                     | T                   | CGT               | CGA             | R              | R            | S               |
| geneB     | contig1  |        924 |       -1 |            82 |              1 | C                     | G                   | GGT               | GCT             | G              | A            | NS              |
| geneB     | contig1  |        913 |       -1 |            93 |              0 | G                     | A                   | CCT               | TCT             | P              | S            | NS              |
| geneB     | contig1  |        909 |       -1 |            97 |              1 | G                     | A                   | CCC               | CTC             | P              | L            | NS              |
| geneB     | contig1  |        905 |       -1 |           101 |              2 | C                     | A                   | CCG               | CCT             | P              | P            | S               |
| geneB     | contig1  |        903 |       -1 |           103 |              1 | G                     | A                   | CCA               | CTA             | P              | L            | NS              |
| geneB     | contig1  |        902 |       -1 |           104 |              2 | T                     | A                   | CCA               | CCT             | P              | P            | S               |


## Codon-Level Events (grouped)

| gene_id   |   codon_start_index |   k | positions   | pos_in_codon   | ancestral_codon   | derived_codon   | ancestral_aa   | derived_aa   | mutation_type   |
|:----------|--------------------:|----:|:------------|:---------------|:------------------|:----------------|:---------------|:-------------|:----------------|
| geneA     |                   9 |   1 | 110         | 2              | GCA               | GCG             | A              | A            | S               |
| geneA     |                  12 |   1 | 113         | 2              | GGT               | GGC             | G              | G            | S               |
| geneA     |                  15 |   1 | 116         | 2              | GGC               | GGT             | G              | G            | S               |
| geneA     |                  21 |   1 | 120         | 0              | GGA               | AGA             | G              | R            | NS              |
| geneA     |                  24 |   1 | 124         | 1              | CCT               | CTT             | P              | L            | NS              |
| geneA     |                  27 |   1 | 126         | 0              | CCC               | GCC             | P              | A            | NS              |
| geneA     |                  51 |   2 | 151,152     | 1,2            | GTC               | GCT             | V              | A            | NS              |
| geneA     |                  60 |   2 | 159,161     | 0,2            | TCT               | GCC             | S              | A            | NS              |
| geneA     |                  72 |   1 | 173         | 2              | CGT               | CGG             | R              | R            | S               |
| geneA     |                  75 |   1 | 174         | 0              | CGC               | GGC             | R              | G            | NS              |
| geneA     |                  81 |   1 | 182         | 2              | CGA               | CGC             | R              | R            | S               |
| geneA     |                  87 |   1 | 187         | 1              | GCC               | GAC             | A              | D            | NS              |
| geneA     |                 105 |   1 | 206         | 2              | GGA               | GGT             | G              | G            | S               |
| geneA     |                 108 |   1 | 207         | 0              | CCT               | ACT             | P              | T            | NS              |
| geneA     |                 123 |   2 | 222,224     | 0,2            | ACC               | TCA             | T              | S            | NS              |
| geneB     |                  12 |   1 | 994         | 0              | CCC               | GCC             | P              | A            | NS              |
| geneB     |                  15 |   1 | 989         | 2              | CCG               | CCC             | P              | P            | S               |
| geneB     |                  27 |   1 | 977         | 2              | ACG               | ACA             | T              | T            | S               |
| geneB     |                  30 |   1 | 974         | 2              | ACA               | ACG             | T              | T            | S               |
| geneB     |                  57 |   2 | 947,949     | 0,2            | CGT               | TGA             | R              | *            | NS              |
| geneB     |                  81 |   1 | 924         | 1              | GGT               | GCT             | G              | A            | NS              |
| geneB     |                  93 |   1 | 913         | 0              | CCT               | TCT             | P              | S            | NS              |
| geneB     |                  96 |   1 | 909         | 1              | CCC               | CTC             | P              | L            | NS              |
| geneB     |                  99 |   1 | 905         | 2              | CCG               | CCT             | P              | P            | S               |
| geneB     |                 102 |   2 | 902,903     | 1,2            | CCA               | CTT             | P              | L            | NS              |


### Notes & Conventions

- Positions in profiles are contig-0-based; Prodigal coordinates are 1-based.

- Reverse-strand genes: profile alleles are forward-strand, but inserted into coding sequence as complements.

- `mutation_type` reflects the final ancestral vs derived codon translation (NCBI table 11).
