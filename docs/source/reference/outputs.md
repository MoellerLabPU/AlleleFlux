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

### Site-Depth Filter Sidecars

**Path:** `allele_freq_cache/{timepoint}/{mag}_{timepoint}_allele_frequency_site_depth_filter.tsv`

Written next to every allele-frequency cache file: what the per-site depth filter (`quality_control.site_depth_filter`) did to each sample's profile before the cache was built. One row per sample.

| Column | Type | Description |
|--------|------|-------------|
| `sample_id`, `mag_id`, `timepoint` | str | The sample and MAG the row describes |
| `mean_coverage` | float | The sample's mean coverage over covered positions (QC's `average_coverage`), the number each position's depth was compared against |
| `mean_coverage_after_filter` | float | The same mean over the positions that were kept; the drop from `mean_coverage` is how much pile the filter removed from this sample |
| `n_sites_covered` | int | Positions with at least one read before filtering |
| `n_removed_low` | int | Positions dropped for depth below `site_depth` (0 with the default floor of 0) |
| `n_removed_high` | int | Positions dropped for depth above `site_ratio` × `mean_coverage` |
| `n_kept` | int | Positions that entered the cache |

The workflow's `site_depth_filter_summary` rule runs `alleleflux-site-filter-summary` once every cache job has finished and writes `site_depth_filter/site_depth_filter_summary.tsv` (every sidecar stacked into one run-level table) and `site_depth_filter/site_depth_filter_summary_by_mag.tsv` next to it, one row per MAG × timepoint with `n_samples`, `mean_of_mean_coverage`, `mean_of_mean_coverage_after_filter`, `sites_removed_high_total`, `sites_covered_total` and `fraction_removed_high` (the share of covered positions the ceiling removed). Pairwise ANI, baseline presence, dN/dS and the plotting tools read this cache rather than raw profiles, so a position removed for a sample is absent for that sample everywhere downstream.

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

:::{warning}
**Experimental.** Outlier gene detection (binomial and Poisson tests, `alleleflux-outliers`) was not used or benchmarked in the AlleleFlux publication. It is off by default (`use_outlier_detection: False`); these files are only produced if you opt in.
:::

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

## Strain Turnover and Baseline Presence Outputs

Written only when the corresponding `analysis.use_*` flag is on. See the [Strain Turnover and Baseline Presence guide](../usage/strain_turnover_analysis.md).

### Pairwise ANI

**Path:** `pairwise_ani/{mag}_pairwise_ani.tsv`, `pairwise_ani/{mag}_pairwise_ani_samples.tsv`

One row per compared sample pair.

| Column | Description |
|--------|-------------|
| `compared_bases_count`, `percent_genome_compared`, `coverage_overlap` | Positions with ≥ `min_cov` reads in both samples, and that as a fraction of the genome |
| `consensus_SNPs`, `conANI` | Positions whose majority base differs, and 1 − that / compared |
| `population_SNPs`, `popANI` | Positions sharing no allele, and 1 − that / compared |
| `min_cov` | The depth gate the metric was computed with |
| `subjectID_1/2`, `group_1/2`, `time_1/2`, `replicate_1/2` | Each sample's metadata |

### Strain Turnover

**Path:** `strain_turnover/{mag}_strain_turnover.tsv`, `strain_turnover/{mag}_turnover_rollup.tsv`

One row per subject × transition.

| Column | Description |
|--------|-------------|
| `transition`, `sample_t1`, `sample_t2` | The earlier and later sample of this subject |
| `compared_bases_count`, `percent_genome_compared`, `conANI`, `popANI` | Copied from the pair |
| `strain_replacement` | popANI < `pop_threshold`; `<NA>` when undetermined |
| `dominant_strain_change` | conANI < `con_threshold`; `<NA>` when undetermined |
| `background` | `stable`, `strain_replacement`, `dominant_strain_change` or `undetermined` |
| `frequency_shift`, `min_compared`, `pop_threshold`, `con_threshold`, `min_cov` | Provenance: the thresholds the verdicts were called with |

### Replacement Classification

**Path:** `strain_turnover/replacement_classification.tsv`

One row per MAG × group × transition × metric (`strain_replacement` and `dominant_strain_change` stacked in a `metric` column).

| Column | Description |
|--------|-------------|
| `n_mice_with_call`, `n_mice_undetermined`, `n_mice_changed` | Mice with a verdict, without one, and with a True verdict |
| `majority_mice_changed`, `any_mouse_changed`, `all_mice_changed`, `no_mouse_changed` | Mouse-block flags. **Blank** when fewer than `min_voters` mice had a verdict: a blank means "no verdict", never "not changed" |
| `strain_status` | **The column to filter on.** Under the default `vote_rule: majority`: `replaced` (at least `min_voters` mice with a verdict, more than half changed), `not_replaced` (at least `min_voters`, half or fewer changed), `too_few_voters` (1 to `min_voters`-1), `no_voters` (none), and `tie` (exactly half, only with `tie: unresolved`). With `vote_rule: any` or `all` the line between the first two moves; the last three do not |
| `n_replicates_with_call`, `n_replicates_changed`, `majority_replicates_changed`, `all_replicates_changed` | Replicate block; whether a replicate counts as changed follows `replicate_rule` (by default: if the mean ANI of its mice with a verdict is below the threshold). Flags blank under the same floor, counted in replicates |
| `strain_status_replicates` | The same status word with replicates as the voters |
| `min_voters`, `vote_rule`, `tie_rule`, `replicate_rule` | The settings that produced the row |

A key is in one of three situations (replaced, checked and not replaced, not
checkable) and a True/False flag can hold only two. Keep
`strain_status == "not_replaced"` to select verified MAGs; testing a flag for
`True` would let every unchecked MAG through.

### Baseline Presence

**Path:** `baseline_presence/{comparison}_{family}_{statistic}_baseline_presence.tsv.gz` (long), `..._summary.tsv`, `..._rising.tsv`

Long table: one row per significant site × allele × sample, with `allele_reads`, `total_reads`, `detection_threshold_reads`, `allele_frequency`, `allele_status` (`present`, `below_detection`, `absent`, `not_covered`; `not_covered` covers both "too few reads" and "position removed for this sample by the site-depth filter", since neither leaves a cache row), `origin_in_own_mouse`, and optionally `strain_background`.

Summary: one row per site × allele. Sample-count columns (`origin_any_mouse`, `n_{earlier}_samples_allele_present`, `n_{earlier}_samples_covered`, `n_replicates_with_allele_at_{earlier}`, `n_mice_standing_variation`, `n_mice_de_novo_candidate`, ...) apply the presence rule; read-count columns (`total_reads_{tp}`, `allele_reads_{tp}`, `allele_frequency_{tp}`) are unfiltered sums over every sample at that timepoint. Timepoint names in column headers are the comparison's own labels.

Rising table: the summary's columns plus `rose_in`, `frequency_change`, `n_replicates`, `n_mice`; one row per (site, allele, group the allele rose in). At a biallelic site both alleles are candidates and both are listed in the summary, but only one rose; the one that fell was the common allele at baseline and always reads standing. The rising table leaves it out, and checks end-point presence in the `rose_in` group only. Count standing vs de novo from this table.

### Rising Origin

**Path:** `baseline_presence/rising_origin/rising_origin_summary.tsv`, `rising_origin_alleles.tsv`

Written by `alleleflux-rising-origin` once every comparison's rising table exists. Each rising allele gets one origin:

| origin | meaning |
|---|---|
| `standing_variation` | present in at least one covered earlier sample |
| `de_novo` | earlier samples covered, the allele in none of them |
| `de_novo_below_detection` | a few earlier reads, under the detection bar |
| `new_strain` | (with strain verdicts) de novo, but in a MAG whose strain was replaced in that comparison x group: the allele most likely arrived with the new strain |
| `not_present_later` | no covered later sample of the group shows it |
| `baseline_not_covered` | no earlier sample deep enough to judge |

The first four are verdicts. With strain turnover on, every allele is joined to its MAG's verdict for that period: a within-group allele reads the group it rose in, and a between-group allele reads both groups and takes the dirtier arm. A MAG never classified is `unverified` (`strain_status = not_in_strain_table`).

Summary: one row per floor × comparison × group × strain class × origin, zero counts included, with `n_alleles`, `n_mags`, `n_rising`, `n_called`, `pct_of_all` (of every rising allele in the block) and `pct_of_called` (of the four verdicts). Strain classes are each class present plus `all`; without strain verdicts only `all`, with no floor. Comparisons and groups are never pooled.

## File Format Notes

- Most files are gzip-compressed TSV (`.tsv.gz`)
- Position numbering is **0-based**
- Missing values: `NaN` or empty string
- p-values: [0, 1] range; significant sites typically p < 0.05
- Allele frequencies: [0, 1] range (proportion of reads)

See also: [Interpreting Results](../usage/interpreting_results.md), [CLI Reference](cli_reference.md)
