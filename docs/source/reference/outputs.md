# Output Files Reference

AlleleFlux generates structured outputs organized by analysis type and data type.

## Output Directory Structure

```text
{root_dir}/
└── {data_type}/                    # "single" or "longitudinal"
    ├── profiles/                   # Sample profiles (per MAG)
    ├── inputMetadata/              # MAG-sample mappings
    ├── QC/                         # Quality control metrics
    ├── eligibility_table_*.tsv     # MAG test eligibility
    ├── allele_analysis/            # Allele frequencies
    ├── significance_tests/         # Statistical test results
    │   ├── two_sample_unpaired_*/
    │   ├── two_sample_paired_*/
    │   ├── single_sample_*/
    │   ├── lmm_*/
    │   └── cmh_*/
    ├── scores/                     # Parallelism & divergence scores
    │   ├── intermediate/MAG_scores_*/
    │   └── processed/
    │       ├── combined/           # MAG-level summaries
    │       └── gene_scores_*/      # Gene-level summaries
    └── outlier_genes/              # Outlier gene detection
```

## Core Output Files

### Profile Files

**Path:** `profiles/{sample}/{sample}_{mag}_profiled.tsv.gz`

Base-level coverage and allele counts per sample-MAG pair.

| Column | Type | Description |
|--------|------|-------------|
| `contig` | str | Contig identifier |
| `position` | int | 0-based position |
| `ref_base` | str | Reference base (A/C/G/T/N) |
| `total_coverage` | int | Total read depth |
| `A`, `C`, `G`, `T` | int | Base counts |
| `gene_id` | str | Overlapping gene (null if intergenic) |

### Quality Control Files

**Path:** `QC/QC_{timepoints}-{groups}/{mag}_qc.tsv`

Sample-level QC metrics.

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | str | Sample identifier |
| `breadth_of_coverage` | float | Fraction of genome covered (0-1) |
| `mean_coverage` | float | Average depth |
| `passed_breadth` | bool | Whether sample passed QC |

### Eligibility Table

**Path:** `eligibility_table_{timepoints}-{groups}.tsv`

Determines which MAGs qualify for each statistical test.

| Column | Description |
|--------|-------------|
| `mag_id` | MAG identifier |
| `unpaired_test_eligible` | Eligible for unpaired tests and LMM |
| `paired_test_eligible` | Eligible for paired tests and CMH |
| `single_sample_eligible_{group}` | Per-group single-sample eligibility |

### Allele Frequency Files

**Path:** `allele_analysis/allele_analysis_{timepoints}-{groups}/{mag}_allele_frequency_*.tsv.gz`

Position-level allele frequencies across samples.

| Column | Type | Description |
|--------|------|-------------|
| `contig` | str | Contig identifier |
| `position` | int | 0-based position |
| `ref_base` | str | Reference base |
| `{sample}_allele_freq` | float | Allele frequency in sample (0-1) |
| `{sample}_alt_allele` | str | Most common non-reference allele |
| `{sample}_coverage` | int | Read depth |
| `gene_id` | str | Overlapping gene |

## Statistical Test Results

### Two-Sample Tests

**Paths:**
- `significance_tests/two_sample_unpaired_{timepoints}-{groups}/{mag}_*.tsv.gz`
- `significance_tests/two_sample_paired_{timepoints}-{groups}/{mag}_*.tsv.gz`

| Column | Type | Description |
|--------|------|-------------|
| `contig`, `position` | str, int | Genomic location |
| `gene_id` | str | Overlapping gene |
| `tTest_p_value` | float | T-test p-value |
| `mannwhitneyu_p_value` | float | Mann-Whitney U p-value |
| `mean_diff` | float | Mean allele frequency difference |
| `cohen_d` | float | Effect size |

### Single-Sample Test

**Path:** `significance_tests/single_sample_{timepoints}-{groups}/{mag}_*.tsv.gz`

Tests deviation from reference within each group.

| Column | Description |
|--------|-------------|
| `avg_allele_freq_{group}` | Mean allele frequency in group |
| `tTest_p_value_{group}` | One-sample t-test p-value |

### CMH Test

**Path:** `significance_tests/cmh_{timepoints}-{groups}/{mag}_*.tsv.gz`

Cochran-Mantel-Haenszel test stratified by replicate/timepoint.

| Column | Description |
|--------|-------------|
| `cmh_p_value` | CMH test p-value |
| `mode` | `across-time` or `across-group` |
| Stratum columns | Allele counts per stratum (replicate or timepoint) |

### LMM Test

**Path:** `significance_tests/lmm_{timepoints}-{groups}/{mag}_*.tsv.gz`

Linear mixed-effects model for longitudinal data.

| Column | Description |
|--------|-------------|
| `lmm_p_value` | Fixed-effect p-value |
| `coefficient` | Estimated effect size |

## Score Files

### MAG-Level Scores

**Path:** `scores/processed/combined/{timepoints}-{groups}/{test_type}_mag_scores.tsv.gz`

Parallelism and divergence scores per MAG.

| Column | Description |
|--------|-------------|
| `mag_id` | MAG identifier |
| `num_significant_sites` | Count of significant positions |
| `parallelism_score` | Score for within-group consistency (0-1) |
| `divergence_score` | Score for between-group difference (0-1) |
| `combined_score` | `parallelism_score × divergence_score` |

### Gene-Level Scores

**Path:** `scores/processed/gene_scores_{timepoints}-{groups}/{test_type}_gene_scores.tsv.gz`

Scores aggregated by gene.

| Column | Description |
|--------|-------------|
| `mag_id`, `gene_id` | Identifiers |
| `gene_parallelism_score` | Mean parallelism across gene positions |
| `gene_divergence_score` | Mean divergence across gene positions |
| `num_significant_positions` | Count of significant sites in gene |

### Outlier Gene Files

**Path:** `outlier_genes/{timepoints}-{groups}/{test_type}_outlier_genes.tsv.gz`

Genes with exceptionally high scores (potential adaptive targets).

| Column | Description |
|--------|-------------|
| `mag_id`, `gene_id` | Identifiers |
| `parallelism_score` | Gene-level parallelism score |
| `outlier_type` | `parallelism`, `divergence`, or `combined` |
| `z_score` | Standard deviations from MAG mean |

## dN/dS Analysis Outputs

Generated by `alleleflux-dnds-from-timepoints` (see {doc}`../usage/dnds_analysis`).

**Codon Events:** `{mag}_codon_events_ng86.tsv.gz` – Path-averaged S/NS counts per codon

**Gene Summary:** `{mag}_gene_summary_ng86.tsv.gz` – dN/dS ratios per gene

**MAG Summary:** `{mag}_mag_summary_ng86.tsv.gz` – Overall MAG dN/dS

**Global Summary:** `{mag}_global_summary_ng86.tsv` – Aggregate statistics

Key columns:
- `dN_dS`: dN/dS ratio (>1 = positive selection, <1 = purifying)
- `potential_S`, `potential_N`: Expected synonymous/non-synonymous sites
- `observed_S`, `observed_N`: Fractional observed counts (path-averaged)
- `k`: Number of positions changed in codon (1, 2, or 3)

## File Format Notes

- Most files are gzip-compressed TSV (`.tsv.gz`)
- Position numbering is **0-based**
- Missing values: `NaN` or empty string
- p-values: [0, 1] range; significant sites typically p < 0.05
- Allele frequencies: [0, 1] range (proportion of reads)

See also: {doc}`../usage/interpreting_results`, {doc}`cli_reference`
