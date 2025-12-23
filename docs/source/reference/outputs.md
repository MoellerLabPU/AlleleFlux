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

:::{seealso}
For detailed information about all statistical tests and score calculation formulas, see [Statistical Tests Reference](statistical_tests.md).
:::

**How Scores Are Calculated:**

For **most tests** (two-sample, LMM, single-sample), the score represents the **percentage of significant sites**:

```text
Score (%) = (Significant Sites / Total Sites) × 100
```

For **CMH tests**, the score uses **differential significance** between timepoints (see CMH section below).

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
| ------ | ----------- |
| `p_value_CMH` | CMH test p-value |
| `time` | Timepoint identifier (for longitudinal data) |
| `mode` | `across-time` or `across-group` |
| Stratum columns | Allele counts per stratum (replicate or timepoint) |

**CMH Score Calculation (Differential Significance):**

Unlike other tests, CMH scores measure **sites that become significant only at the focus timepoint**:

1. **Identify common sites**: Positions present in results from BOTH timepoints
2. **Find differential sites**: Sites where `p_value_CMH < threshold` at focus timepoint BUT NOT at the other timepoint
3. **Calculate percentage**: `Score (%) = (Differential Sites / Common Sites) × 100`

**Mathematical formula:**

```text
Let S_focus = sites significant at focus timepoint
Let S_other = sites significant at other timepoint
Let S_common = sites in both timepoints

Differential = (S_focus - S_other) ∩ S_common
Score = |Differential| / |S_common| × 100
```

**Interpretation:**

- High scores indicate strong timepoint-specific selection at the focus timepoint
- The focus timepoint (typically the later/derived state) determines directionality
- Example: If `focus: post`, the score measures sites that became significant from pre→post

### LMM Test

**Path:** `significance_tests/lmm_{timepoints}-{groups}/{mag}_*.tsv.gz`

Linear mixed-effects model for longitudinal data.

| Column | Description |
|--------|-------------|
| `lmm_p_value` | Fixed-effect p-value |
| `coefficient` | Estimated effect size |

## Score Files

### MAG-Level Scores

**Path:** `scores/processed/combined/MAG/scores_{test_type}-{timepoints}-{groups}*.tsv`

Evolutionary significance scores per MAG for each statistical test.

**Standard Tests (two-sample, LMM, single-sample, etc.):**

| Column | Description |
| ------ | ----------- |
| `MAG_ID` | MAG identifier |
| `total_sites_per_group_{test}` | Total genomic positions analyzed |
| `significant_sites_per_group_{test}` | Number of sites with p < threshold |
| `score_{test} (%)` | Percentage of significant sites: `(significant/total) × 100` |
| Taxonomy columns | Domain, phylum, class, order, family, genus, species |
| `grouped_by` | Grouping level (e.g., "MAG_ID") |

**CMH Scores (special calculation):**

| Column | Description |
| ------ | ----------- |
| `MAG_ID` | MAG identifier |
| `focus_timepoint` | Which timepoint was designated as focus (derived state) |
| `total_sites_per_group_CMH` | Number of common sites across both timepoints |
| `significant_sites_per_group_CMH` | Differential significant sites (focus only) |
| `score_CMH (%)` | Percentage: `(differential sites / common sites) × 100` |
| Taxonomy columns | Domain, phylum, class, order, family, genus, species |
| `grouped_by` | Grouping level |

**Understanding the scores:**

- **Higher scores** indicate more genomic positions showing evolutionary signatures
- **Standard tests**: Direct measure of selection strength across the genome
- **CMH scores**: Measure of timepoint-specific directional changes

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

Generated by `alleleflux-dnds-from-timepoints` (see [dN/dS Analysis Guide](../usage/dnds_analysis.md)).

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

See also: [Interpreting Results](../usage/interpreting_results.md), [CLI Reference](cli_reference.md)
