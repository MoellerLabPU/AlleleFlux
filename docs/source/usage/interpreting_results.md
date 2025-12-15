# Interpreting Results

This guide explains how to interpret the results produced by AlleleFlux.

## Output Structure

AlleleFlux organizes results by analysis type:

```text
output/
├── profiles/               # Per-sample allele counts
├── metadata/               # Per-MAG sample metadata
├── QC/                     # Quality control results
├── eligibility_table_*.tsv # MAG eligibility for tests
├── allele_analysis/        # Allele frequency analysis
├── significance_tests/     # Statistical test results
│   ├── two_sample_unpaired/
│   ├── two_sample_paired/
│   ├── single_sample/
│   ├── lmm/
│   └── cmh/
├── scores/
│   ├── intermediate/        # Per-MAG scores
│   └── processed/
│       ├── combined/        # Aggregated MAG/taxa scores
│       └── gene_scores/     # Gene-level scores
└── outliers/               # High-scoring outlier genes
```

## Key Files

**1. Eligibility Table**

`eligibility_table_{timepoints}-{groups}.tsv` - Which MAGs qualify for each test based on coverage/samples

**2. Statistical Tests** (in `significance_tests/`)

Per-MAG files with p-values and test statistics:
\- `{mag}_two_sample_unpaired.tsv.gz` - Unpaired group comparisons
\- `{mag}_lmm.tsv.gz` - Linear mixed models
\- `{mag}_cmh.tsv.gz` - Cochran-Mantel-Haenszel tests

Key columns: `contig`, `position`, `gene_id`, `p_value_{test}`, `q_value_{test}`

**3. Scores** (in `scores/processed/combined/`)

- `scores_{test}-{tp}-{gr}-MAGs.tsv` - MAG-level parallelism/divergence scores
- `scores_{test}-{tp}-{gr}-{taxon}.tsv` - Taxonomic aggregations (phylum to species)

**4. Gene Scores** (in `scores/processed/gene_scores/`)

- `{mag}_{test}_gene_scores_individual.tsv` - Per-gene scores
- `{mag}_{test}_outlier_genes.tsv` - High-scoring genes under selection

## Score Interpretation

**Parallelism Score** (0-100%)
Measures consistent allele changes across replicates within a group. High scores → deterministic evolution (not random drift).

**Divergence Score** (0-100%)
Quantifies allele frequency differences between groups. High scores → differential selection between conditions.

**CMH Test**
Detects parallel allele changes across timepoints while controlling for individual variation. Particularly powerful for longitudinal studies.

## File Format Details

**Profile files** (`profiles/{sample}_{mag}_profiled.tsv.gz`):
`contig`, `position`, `ref_base`, `total_coverage`, `A`, `C`, `G`, `T`, `gene_id`

**Statistical test results** (`significance_tests/{test}/{mag}_{test}.tsv.gz`):
`contig`, `position`, `gene_id`, `p_value_{test}`, `q_value_{test}`

**Gene scores** (`scores/processed/gene_scores/{mag}_{test}_gene_scores_individual.tsv`):
`gene_id`, `total_sites`, `significant_sites`, `score_%`

**Outliers** (`outliers/{mag}_{test}_outlier_genes.tsv`):
`gene_id`, `gene_score_%`, `mag_score_%`, `p_value_binomial`, `p_value_poisson`

## Analysis Workflow

**Step 1: Check Eligibility**

```bash
cat eligibility_table_pre_post-treatment_control.tsv
```

Identify MAGs with sufficient coverage for statistical tests.

**Step 2: Examine Scores**

```bash
# MAG-level scores
head scores_two_sample_unpaired-pre_post-treatment_control-MAGs.tsv

# Taxonomic aggregation (family level)
head scores_two_sample_unpaired-pre_post-treatment_control-family.tsv
```

Focus on MAGs/taxa with high parallelism or divergence scores.

**Step 3: Investigate Genes**

```bash
# Gene scores for a high-scoring MAG
head MAG123_two_sample_unpaired_gene_scores_individual.tsv

# Outlier genes
head MAG123_two_sample_unpaired_outlier_genes.tsv
```

Identify candidate genes under strong selection.

**Step 4: Compare Tests**

Check consistency across statistical approaches (two-sample, LMM, CMH). Genes significant in multiple tests are most robust.

**Step 5: Functional Analysis**

- Annotate outlier genes (KEGG, COG, Pfam)
- Check biological relevance to experimental conditions
- Consider genomic context (operons, mobile elements)

## Troubleshooting

**No results / empty files**
\- Check eligibility table: MAGs may not meet `min_sample_num` or `breadth_threshold`
\- Verify input file paths in configuration
\- Check log files in `logs/` directory

**Low scores across all MAGs**
\- Insufficient selective pressure or inappropriate timepoints
\- Try lowering `p_value_threshold` (e.g., 0.1 instead of 0.05)
\- Check if experimental conditions are strong enough

**Inconsistent results between tests**
\- LMM is sensitive to experimental design complexity
\- Two-sample tests affected by unbalanced groups
\- CMH best for detecting consistent directional changes
\- Use multiple tests for robust conclusions

**Missing gene IDs**
\- Ensure Prodigal predictions match reference FASTA contig names
\- Verify `prodigal_path` in configuration
\- Check gene FASTA headers match contig naming

For visualization of results, see {doc}`visualization_guide`.
