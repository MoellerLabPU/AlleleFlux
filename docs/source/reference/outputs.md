# Output Files Reference

AlleleFlux generates structured outputs organized by analysis type and data type.

## Output Directory Structure

```text
{root_dir}/
└── {data_type}/                    # "single" or "longitudinal"
    ├── profiles/                   # Sample profiles (per MAG)
    ├── inputMetadata/              # MAG-sample mappings, one folder per timepoint
    │   └── inputMetadata_{timepoint}/
    ├── QC/                         # Quality control metrics, one folder per timepoint
    │   └── QC_{timepoint}/
    ├── eligibility_table_*.tsv     # MAG test eligibility, one per comparison
    ├── allele_freq_cache/          # Per-sample allele frequencies, one folder per timepoint
    │   └── {timepoint}/
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

### Sample Metadata Files

**Path:** `inputMetadata/inputMetadata_{timepoint}/{mag}_metadata.tsv`, plus `summary_mag_index.tsv` in the same folder

Built once per **single timepoint**, not per comparison and not per group pair. With comparisons `pre_post` and `pre_end` there are three folders (`pre`, `post`, `end`), and the `pre` samples are listed once. For single-timepoint data the comparison label already is one timepoint, so there is one folder per comparison label.

Each file lists the profiled samples of one MAG at that timepoint, for every group named in `groups_combinations` (plus the sheet's other groups when they are read as baseline witnesses, see `cohort_witnesses` in the [Strain Turnover and Baseline Presence guide](../usage/strain_turnover_analysis.md)): `sample_id`, `subjectID`, `time` (longitudinal data only), `group`, `replicate`, `file_path` (the profile), `bam_path`.

Every MAG that has a profile somewhere in the run gets a file in every timepoint folder. A MAG with no profiled sample at a timepoint gets a file with the header and no rows. So "no sample here" is a file you can open, and a missing file means a step failed.

`summary_mag_index.tsv` stacks the rows of every MAG of the folder into one table (`mag_id` first).

### Quality Control Files

**Path:** `QC/QC_{timepoint}/{mag}_QC.tsv`

Built once per **single timepoint**, from that timepoint's metadata folder. One row per sample listed in that metadata folder, passing or not, all groups together. The table holds facts about one sample and one MAG only; nothing in it depends on a comparison or a group pair.

| Column | Type | Description |
|--------|------|-------------|
| `sample_id`, `MAG_ID`, `file_path` | str | The sample, the MAG and the profile that was read |
| `group`, `subjectID`, `replicate` | str | Copied from the sample sheet |
| `genome_size` | int | Length of the MAG in bases |
| `breadth` | float | Fraction of the genome with at least one read (0-1) |
| `breadth_threshold_passed`, `breadth_fail_reason` | bool, str | The breadth check against `quality_control.breadth_threshold` |
| `average_coverage` | float | Mean depth over the covered positions only |
| `average_coverage_including_zeros` | float | Summed depth / genome size; the number `quality_control.coverage_threshold` tests |
| `coverage_threshold_passed`, `coverage_fail_reason` | bool, str | **The sample's QC verdict.** True only when the breadth check passed too |
| `median_coverage`, `median_coverage_including_zeros`, `coverage_std`, `coverage_std_including_zeros`, `length_weighted_coverage` | float | Further depth statistics |
| `time` | str | The timepoint (longitudinal data only) |

A MAG with no profiled sample at the timepoint gets a table with the header and no rows.

Three summary files sit next to the per-MAG tables and cover the same single timepoint: `ALL_MAGs_QC_group_summary.tsv`, `ALL_MAGs_QC_group_time_summary.tsv` (longitudinal data only) and `ALL_MAGs_QC_overall_summary.tsv`.

Which subjects were readable at both timepoints of a comparison, and how many replicates each group then has, are **not** in the QC tables. Those are questions about a comparison, and the eligibility step answers them (next section).

### Eligibility Table

**Path:** `eligibility_table_{timepoints}-{groups}.tsv`

One table per comparison (`{timepoints}` is `{earlier}_{later}` for longitudinal data, the timepoint itself for single data) and group pair. It says which MAGs qualify for each statistical test.

| Column | Description |
|--------|-------------|
| `MAG_ID` | MAG identifier |
| `replicates_per_group_{group}` | Replicates of that group among the samples that count (see below) |
| `paired_replicates_per_group_{group}` | Of those, the replicates that both tested groups have |
| `unpaired_test_eligible` | Both groups have at least `min_sample_num` replicates: unpaired tests and LMM |
| `paired_test_eligible` | Both groups have at least `min_sample_num` paired replicates: paired tests and CMH |
| `single_sample_eligible_{group}` | That group has at least `min_sample_num` replicates (longitudinal data; empty for single data) |

**How the counts are made.** The eligibility step reads the QC folders of the comparison's timepoints (two for longitudinal data, one for single) and works out the pairing itself:

1. A sample **counts** when it passed QC (`coverage_threshold_passed`) and, for longitudinal data, its subject has a passing sample at **both** timepoints of the comparison.
2. Only the two tested groups are kept.
3. Replicates are counted per group among the samples that count. A replicate is *paired* when both tested groups have it.

Small example, longitudinal, comparison `pre_post`, groups `A` and `B`, `min_sample_num` 2, replicate = cage:

| subject | group | cage | passed QC at `pre` | passed QC at `post` | counts? |
|---|---|---|---|---|---|
| a1 | A | c1 | yes | yes | yes |
| a2 | A | c2 | yes | yes | yes |
| b1 | B | c3 | yes | yes | yes |
| b2 | B | c4 | yes | no | no (not readable at `post`) |

| step | result |
|---|---|
| replicates per group | A: 2 cages (c1, c2); B: 1 cage (c3) |
| paired replicates | 0 (no cage holds both groups) |
| unpaired test | not eligible: B is below the floor of 2 |
| paired test | not eligible: 0 paired replicates |
| single-sample test | A eligible, B not |

A MAG appears in the table only when both tested groups have at least one sample that counts; a MAG with samples at only one of the two timepoints is left out. Rows are in MAG name order.

### Allele-Frequency Cache

**Path:** `allele_freq_cache/{timepoint}/{mag}_{timepoint}_allele_frequency.parquet`

One Parquet file per MAG per **single timepoint**, holding one row per (contig, position, sample) for every QC-passing sample of that timepoint, all groups together. It is built from the MAG's QC table of the same timepoint and is read by every comparison that names the timepoint. See [Allele Frequency Cache Architecture](../usage/allele_freq_cache_architecture.md).

Which (MAG, timepoint) pairs get a cache:

| Run | Caches built |
|---|---|
| Any run | Each MAG that is eligible for at least one test in a comparison, at that comparison's timepoints |
| `baseline_scope: cohort`, longitudinal data, every sample dated, baseline presence or gene de novo switched on | Each such MAG at **every witness timepoint**: every compared timepoint, plus (under `cohort_witnesses: sheet`, the default) every other timepoint of the sample sheet |

The second row exists because, under the cohort scope, a sample at *any* timepoint can show that an allele was already there (see [Which samples are the baseline](../usage/strain_turnover_analysis.md#which-samples-are-the-baseline-baseline_scope)). Example: comparisons `pre_post` and `pre_end`, the sheet also holds timepoint `mid`, and MAG `M2` is eligible in `pre_post` only. Without the cohort condition `M2` has caches at `pre` and `post`. With it, `M2` has caches at `pre`, `post`, `end` and `mid`.

**Zero-row caches.** A MAG can have no QC-passing sample at one of those timepoints. Its cache is then a Parquet file with the normal columns and **no rows**, and its sidecar (see [Site-Depth Filter Sidecars](#site-depth-filter-sidecars) below) has the header only. Baseline presence and the gene-level de novo analysis read such a file as "no sample of this MAG at this timepoint". A cache file that is missing altogether stops them with an error naming the path.

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

Written next to every allele-frequency cache file: what the per-site depth filter (`quality_control.site_depth_filter`) did to each sample's profile before the cache was built. One row per sample; header only next to a zero-row cache.

| Column | Type | Description |
|--------|------|-------------|
| `sample_id`, `mag_id`, `timepoint` | str | The sample and MAG the row describes |
| `mean_coverage` | float | The sample's mean coverage over covered positions (QC's `average_coverage`), the number each position's depth was compared against |
| `mean_coverage_after_filter` | float | The same mean over the positions that were kept; the drop from `mean_coverage` is how much pile the filter removed from this sample |
| `n_sites_covered` | int | Positions with at least one read before filtering |
| `n_removed_low` | int | Positions dropped for depth below `site_depth` (0 with the default floor of 0) |
| `n_removed_high` | int | Positions dropped for depth above `site_ratio` × `mean_coverage` |
| `n_kept` | int | Positions that entered the cache |

The workflow's `site_depth_filter_summary` rule runs `alleleflux-site-filter-summary` once every cache job has finished and writes `site_depth_filter/site_depth_filter_summary.tsv` (every sidecar stacked into one run-level table) and `site_depth_filter/site_depth_filter_summary_by_mag.tsv` next to it, one row per MAG × timepoint with `n_samples`, `mean_of_mean_coverage`, `mean_of_mean_coverage_after_filter`, `sites_removed_high_total`, `sites_covered_total` and `fraction_removed_high` (the share of covered positions the ceiling removed). A header-only sidecar holds no sample, so a zero-row cache adds no row to either table. Pairwise ANI, baseline presence, dN/dS and the plotting tools read this cache rather than raw profiles, so a position removed for a sample is absent for that sample everywhere downstream.

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
| `num_pairs` | Replicates the test used at this site (see below) |
| `time` | Timepoint identifier (for longitudinal data) |
| `mode` | `across-time` or `across-group` |
| Stratum columns | Allele counts per stratum (replicate or timepoint) |

**Which replicates count at a site.** A replicate is a pair at a site only when **both** groups have a row for it at that site (for the across-time test: both timepoints). "A row" means a sample of that replicate has the position in its allele-frequency data, that is, it has reads there. A replicate with reads in one group and no row in the other is not a pair. The site is tested only when at least `min_sample_num` pairs remain. Example with replicates 1 to 4 and `min_sample_num` 4: replicate 4 has reads in the treatment group and no row in the control group, so the site has 3 pairs and is not tested. See [Statistical Tests Reference](statistical_tests.md).

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
- Example: For the pair `earlier: pre`, `later: post`, the score measures sites that became significant from pre→post

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

Long table: one row per significant site × allele × sample, with `allele_reads`, `total_reads`, `detection_threshold_reads`, `allele_frequency`, `allele_status` (`present`, `below_detection`, `absent`, `not_covered`; `not_covered` covers both "too few reads" and "position removed for this sample by the site-depth filter", since neither leaves a cache row), `origin_in_own_mouse`, and optionally `strain_background`. Two sample-context columns, `collection_date` (ISO `YYYY-MM-DD`, e.g. `2017-01-10`) and `cohort`, are copied from the metadata when it has them and are blank otherwise (see [inputs](inputs.md), optional columns).

Summary: one row per site × allele. Sample-count columns (`origin_any_mouse`, `n_{earlier}_samples_allele_present`, `n_{earlier}_samples_covered`, `n_replicates_with_allele_at_{earlier}`, `{earlier}_mice_allele_present`, `{later}_mice_allele_present` (the carriers), `n_mice_standing_variation`, `n_mice_de_novo_candidate`, ...) apply the presence rule; read-count columns (`total_reads_{tp}`, `allele_reads_{tp}`, `allele_frequency_{tp}`) are unfiltered sums over every sample at that timepoint. Two reporting-only columns say how hard the baseline was looked at (no verdict depends on them): `n_baseline_reads`, the reads at the position over the *covered* baseline samples only, and `n_carriers_own_baseline_absent`, the carriers whose own earlier sample was covered and clean, the strictest per-mouse evidence of a new allele. Timepoint names in column headers are the comparison's own labels.

Five witness-count columns follow `n_carriers_own_baseline_absent` in the summary and in the rising table. They are **always in the header**. They are filled only when the run names its pre-treatment timepoints (`analysis.baseline_presence.pre_treatment_timepoints`, dated `cohort` scope only) and are empty otherwise. A *witness* is a covered baseline sample of the row; it is *pre-treatment* when its timepoint label is in the list and *on-treatment* when it is not.

| Column | Counts |
|---|---|
| `n_pre_treatment_witnesses_covered` | Pre-treatment witnesses |
| `n_pre_treatment_witnesses_allele_present` | Of those, the ones carrying the allele (presence rule) |
| `n_pre_treatment_witnesses_below_detection` | Of those, the ones with a few reads of the allele under the detection bar. They do not carry it, so they never make it standing; a `first_seen_after_treatment_began` label with a count here is a weaker call |
| `n_on_treatment_witnesses_covered` | On-treatment witnesses |
| `n_on_treatment_witnesses_allele_present` | Of those, the ones carrying the allele |

The two `covered` counts add up to `n_{earlier}_samples_covered` and the two `allele_present` counts to `n_{earlier}_samples_allele_present`. With the list set, `origin_any_mouse` can also read `first_seen_after_treatment_began` (see [Three labels](../usage/strain_turnover_analysis.md#three-labels-pre_treatment_timepoints)).

Rising table: the summary's columns plus `rose_in`, `frequency_change`, `n_replicates`, `n_mice`; one row per (site, allele, group the allele rose in). At a biallelic site both alleles are candidates and both are listed in the summary, but only one rose; the one that fell was the common allele at baseline and always reads standing. The rising table leaves it out, and checks end-point presence in the `rose_in` group only. Count standing vs de novo from this table.

### Rising Origin

**Path:** `baseline_presence/rising_origin/rising_origin_summary.tsv`, `rising_origin_alleles.tsv`

Written by `alleleflux-rising-origin` once every comparison's rising table exists. Each rising allele gets one origin:

| origin | meaning |
|---|---|
| `standing_variation` | present in at least one covered earlier sample; in a run that names its pre-treatment timepoints, in at least one covered **pre-treatment** witness |
| `first_seen_after_treatment_began` | only in a run that names its pre-treatment timepoints: present in at least one covered witness, but in none sampled before treatment began. A label, not a proof: the allele may have arisen under treatment, or it may have been there all along in no sampled pre-treatment witness. Never relabelled by the cluster rule or the strain gate |
| `de_novo` | earlier samples covered, the allele in none of them |
| `de_novo_below_detection` | a few earlier reads, under the detection bar |
| `new_strain` | de novo, but the strain gate says the allele most likely arrived with a new strain: the group's across-cage vote says replaced (gate `group`), or more than half of the allele's own carriers changed strain (gate `carrier`) |
| `strain_undetermined` | de novo, but the strain gate could not say: fewer than `min_determined_carriers` carriers with a per-mouse verdict (gate `carrier`), or a group vote that is neither replaced nor not replaced: too few voters, no voters, not in the table, a tie (gate `group`). Never counted as de novo |
| `block_event` | de novo (either flavour) with more than `max_de_novo_neighbours` other de novo alleles within `cluster_window_bp` on the same contig, in the same group, at any age: one imported or recombined chunk, never counted as independent mutations. `block_id` (`{mag}:{contig}:{group}:{k}`) names the block and `n_de_novo_neighbours` the count; both columns are in the allele table for every row |
| `not_present_later` | no covered later sample of the group shows it |
| `baseline_not_covered` | no earlier sample deep enough to judge |

The first five are verdicts and are counted in `n_called`; `block_event`, `strain_undetermined` and the last two are not, so a block or an undetermined strain drops out of `n_called` and lowers the de novo share. An allele labelled `first_seen_after_treatment_began` would read `standing_variation` without the pre-treatment list, so `n_called` is the same number with and without it. Two different letters rising at the *same* position (possible only at a site where three or four bases were tested) are two alleles 0 bp apart and form a block: two independent mutations at one base in one group is about as unlikely as it gets, while reads from a related strain mapping with a mismatch there, or an imported chunk carrying a different letter, produce exactly this picture; neither is two clean de novo events. The cluster rule runs before any strain gate, on every de novo candidate, because adjacency is a fact about positions and a gate that ran first could hide a block member's neighbours. With strain turnover on, every allele is joined to its MAG's verdict for that period: a within-group allele reads the group it rose in, and a between-group allele reads both groups and takes the dirtier arm. A MAG never classified is `unverified` (`strain_status = not_in_strain_table`).

**The strain gate** (`--strain_gate`, config `analysis.baseline_presence.strain_gate`) decides which verdict relabels a de novo allele; only the two de novo flavours are ever relabelled. Two verdicts are available and **both are written to every row** whenever their input was given, whatever the gate: under `group` the carrier columns are reporting only, under `carrier` the group columns are. The gate decides which is *read*, so the two can always be compared on the same run.

- the *group* verdict (`strain_status`, `strain_class`): the across-cage vote of the replacement-classification table for the allele's MAG, group and period;
- the *carrier* verdict (`carrier_strain_status`, `n_carriers`, `n_carriers_determined`, `n_carriers_strain_changed`): the mice carrying the allele at the later timepoint (`{later}_mice_allele_present` in the rising table) looked up in the per-mouse `{mag}_strain_turnover.tsv` for the allele's period; a carrier is *determined* when its `dominant_strain_change` verdict is not blank and *changed* when it is true; more than half changed → `new_strain_carriers`, fewer determined than `min_determined_carriers` (default 3) → `strain_undetermined`, otherwise `same_strain_carriers`. A MAG with no turnover table is `strain_undetermined` with a warning.

| gate | relabels from | replaced / new strain carriers | not replaced / same strain | anything else |
|---|---|---|---|---|
| `carrier` (default) | the carrier verdict | `new_strain` | kept | `strain_undetermined` |
| `group` | the group verdict | `new_strain` | kept | `strain_undetermined` |
| `both` | both: one undetermined verdict → `strain_undetermined` (removed, not explained); else one replaced / new-strain verdict → `new_strain`; else kept | | | |
| `none` | nothing | | | |

Allele table (`rising_origin_alleles.tsv`): one row per rising allele with its `origin`, the five witness-count columns of its rising-table row (`n_pre_treatment_witnesses_covered` ... `n_on_treatment_witnesses_allele_present`; empty when the run names no pre-treatment timepoints), `block_id`, `n_de_novo_neighbours` and the strain columns described above.

Summary (`rising_origin_summary.tsv`): one row per floor × comparison × group × strain class × origin, zero counts included, with `n_alleles`, `n_mags`, `n_rising`, `n_called`, `pct_of_all` (of every rising allele in the block) and `pct_of_called` (of the verdicts). The `first_seen_after_treatment_began` rows are there only when the run made the split; the command reads that from the rising tables themselves (their witness-count columns are filled). Strain classes are each class present plus `all`; without strain verdicts only `all`, with no floor. Comparisons and groups are never pooled.

## Gene-Level De Novo Accumulation

**Path:** `gene_de_novo/{comparison}/gene_de_novo_tests_{comparison}.tsv` (one per comparison) plus three files per MAG: `{mag}_gene_de_novo_events.tsv`, `{mag}_gene_de_novo_genes.tsv`, `{mag}_gene_de_novo_cages.tsv`

Written by `alleleflux-gene-de-novo` when `use_gene_de_novo` is on (longitudinal only). The cage (`group` × `replicate`) is the replicate throughout; see the [Gene De Novo guide](../usage/gene_de_novo_analysis.md) for the method.

**Events** (`{mag}_gene_de_novo_events.tsv`): one row per candidate rise, a (cage, contig, position, letter) whose cage-average share went up in the mean-change table, with the columns of the rise (`contig`, `position`, `gene_id`, blank between genes, `replicate`, `group`, `allele`, `diff_mean`, `subjectID_count`) and the verdict: `origin` (`de_novo`, `de_novo_below_detection`, `standing_variation`, `first_seen_after_treatment_began`, `block_event`, `new_strain`, `strain_undetermined`, `not_present_later`, `baseline_not_covered`; the same words as [rising origin](#rising-origin), judged per cage; `first_seen_after_treatment_began` appears only when the run names its pre-treatment timepoints), `block_id` (`{mag}:{contig}:{group}|{replicate}:{k}`: blocks never span two cages), `n_de_novo_neighbours`, `carriers` (the cage's later-timepoint mice carrying the letter), `n_carriers`, `carrier_strain_status`, `n_carriers_determined`, `n_carriers_strain_changed`, `n_baseline_samples_covered`, `n_baseline_samples_present`. The two de novo flavours count as events; a `block_event` counts only under `block_rule: once`. `first_seen_after_treatment_began` is a reported label and never an event, so the genes, cages and tests tables are the same with and without the pre-treatment list; the events table gets no extra column for it.

**Genes** (`{mag}_gene_de_novo_genes.tsv`): one row per gene × cage, zero rows kept: `n_de_novo` (counted events in the gene), `n_block_excluded` (block rows dropped by `block_rule: exclude`), `covered_bp` (readable bases of the gene in that cage: the denominator, never the gene length), `gene_length`, `covered_fraction`, `entered` (readable above `min_gene_covered_fraction`; a cage that does not enter reports `covered_bp` 0 and `n_de_novo` 0 so the tests never see it), `rate` (`n_de_novo / covered_bp`, blank when not entered). Overlapping genes: a position inside two genes carries both ids joined by a comma in the events table (`gene_10,gene_9`), and an event there counts for each of the two genes, because each gene's `covered_bp` includes the shared bases. The cage total counts it once, so the genes' counts can sum to more than `n_de_novo_total`.

**Cages** (`{mag}_gene_de_novo_cages.tsv`): one row per cage: `n_de_novo_total` (every counted event of the cage, between genes included) and `covered_bp_total` (readable bases of the whole genome). These are the "fair share" numbers of Test A.

**Tests** (`gene_de_novo_tests_{comparison}.tsv`): one row per MAG × gene over every MAG of the comparison. Per group `{g}`: `n_cages_{g}` (entered cages), `n_de_novo_{g}`, `expected_{g}` (the sum over cages of the cage's genome-wide events × the gene's share of the cage's readable genome), `rate_ratio_vs_genome_{g}`, `dispersion_{g}` (quasi-Poisson scale, floored at 1), `p_value_vs_genome_{g}` (Test A: chi-square with 1 df on the dispersion-scaled likelihood-ratio statistic), `q_value_vs_genome_{g}`, `p_value_raw_vs_genome_{g}` and `q_value_raw_vs_genome_{g}` (the same statistic **without** the dispersion discount, chi-square with 1 df: the textbook Poisson test that treats every cage as an independent draw), then the uncorrected diagnostics `df_resid_vs_genome_{g}` (cages used − 1), `lrt_vs_genome_{g}` (the likelihood-ratio statistic before the discount) and the statistic on **both** reference curves, `p_value_chi2_vs_genome_{g}` / `p_value_f_vs_genome_{g}` (discounted statistic on chi-square with 1 df and on F(1, `df_resid_vs_genome_{g}`); the chi-square one is `p_value_vs_genome_{g}`) and `p_value_raw_chi2_vs_genome_{g}` / `p_value_raw_f_vs_genome_{g}` (undiscounted statistic on both curves; the chi-square one is `p_value_raw_vs_genome_{g}`, the F one has no statistical justification and is a diagnostic only), `note_{g}`. Between the groups: `rate_{a}`, `rate_{b}` (pooled events per readable base), `rate_ratio_{a}_vs_{b}`, `dispersion_between`, `df_resid_between` (cages − 2), `lrt_between`, `p_value_between` (Test B: F(1, `df_resid_between`) on the scaled statistic), `q_value_between`, `p_value_raw_between` and `q_value_raw_between` (undiscounted, chi-square with 1 df; reported even with one cage per group, where the discounted test is `no_residual_df`), the same four both-curve diagnostics as above (`p_value_chi2_between` / `p_value_f_between`, `p_value_raw_chi2_between` / `p_value_raw_f_between`; here the F one is `p_value_between`), `note_between`. The F twins are NaN whenever there is no residual df. A gene whose raw p is small but whose discounted p is not is one whose evidence comes from a single cage. Also `gene_length` and `min_covered_fraction` (the smallest covered fraction among the entered cages).

| note | meaning |
|---|---|
| (blank) | tested |
| `below_cage_floor` | fewer entered cages than `quality_control.min_sample_num` in that group (Test B: in either group); nothing tested |
| `no_events` | the group's cages have no event in the gene but a non-zero expectation; rate ratio 0, no p-value |
| `no_exposure` | the expectation is zero (the group's cages have no event anywhere in the genome, or no readable base in the gene); nothing to compare against |
| `no_residual_df` | Test B with one cage per group: the dispersion cannot be estimated; the rate ratio is reported, no p-value |

Benjamini-Hochberg runs within each of six families (Test A per group, Test B, and the raw twin of each) over the tested genes only, per MAG when `statistics.fdr_group_by_mag_id` is true, otherwise over every gene of the comparison. Untested genes keep a blank q-value; they are never given p = 1.

## Figures and run report

Written by the local rules `plot_scores` and `run_report` when `analysis.figures.enabled` is on (the default), for every comparison that has combined MAG score tables; nothing when `analysis.allele_analysis_only` is on. See the [Score Figures and Run Report guide](../usage/figures_and_report.md).

**Path:** `figures/{comparison}/`, `figures/run_report.html`

| File | Description |
|------|-------------|
| `figures/{comparison}/score_heatmap_{family}.png` | Rows = every MAG tested in any cell of the family (taxonomy order, or tree tip order with `tree_path`), a phylum strip, then one tile per cell (parallelism per group, divergence) coloured by score 0-100 %, with the smallest BH q as stars (`***` < 0.001, `**` < 0.01, `*` < threshold); grey = not tested in that cell |
| `figures/{comparison}/score_boxplot_{family}.png` | One box per cell, one dot per MAG, a grey line joining a MAG's two parallelism scores, `n=` under each box |
| `figures/{comparison}/score_boxplot_{family}_by_phylum.png` | The same, one panel per phylum with at least `min_mags_per_phylum` MAGs; smaller phyla pooled as `other` |
| `figures/{comparison}/score_figures.done` | The figure files written for the comparison (the workflow's sentinel); `.pdf` files alongside or instead of `.png` under `format: both` / `pdf` |
| `figures/run_report.html` | One self-contained page: run and QC settings, the sample-and-MAG funnel per comparison, MAGs tested and MAGs with a significant site (`min_q_value < threshold`) per test cell, the site-depth filter per timepoint, and every PNG figure embedded with a caption |

`{family}` is one of `ttest`, `ttest_unpaired`, `wilcoxon`, `mannwhitney`, `lmm`, `cmh`; a family with no score table is skipped.

## File Format Notes

- Most files are gzip-compressed TSV (`.tsv.gz`)
- Position numbering is **0-based**
- Missing values: `NaN` or empty string
- p-values: [0, 1] range; significant sites typically p < 0.05
- Allele frequencies: [0, 1] range (proportion of reads)

See also: [Interpreting Results](../usage/interpreting_results.md), [CLI Reference](cli_reference.md)
