# Configuration Reference

Complete reference for AlleleFlux configuration options.

## Quick Start

Copy the template configuration:

```bash
cp alleleflux/smk_workflow/config.template.yml my_config.yml
```

Edit `my_config.yml` with your file paths and analysis parameters.

## Core Parameters

**run_name** (optional)

Unique identifier for this analysis.

```yaml
run_name: my_study_2024
```

---

### input

Paths to required input files.

| Parameter | Description |
|-----------|-------------|
| `fasta_path` | Path to the combined reference FASTA file containing all MAG contigs. Header format should be `<MAG_ID>.fa_<contig_ID>`. |
| `prodigal_path` | Path to Prodigal gene predictions (nucleotide FASTA). Gene IDs must match contig IDs in the reference FASTA. |
| `metadata_path` | Path to sample metadata TSV file. Must contain columns: `sample_id`, `bam_path`, `subjectID`, `group`, `replicate`. For longitudinal data, also include `time`. |
| `gtdb_path` | Path to GTDB-Tk taxonomy file (`gtdbtk.bac120.summary.tsv`). Used for taxonomic aggregation of scores. |
| `mag_mapping_path` | Path to contig-to-MAG mapping file (TSV with `contig_name` and `mag_id` columns). |

**Example:**

```yaml
input:
  fasta_path: /path/to/combined_mags.fasta
  prodigal_path: /path/to/prodigal_genes.fna
  metadata_path: /path/to/sample_metadata.tsv
  gtdb_path: /path/to/gtdbtk.bac120.summary.tsv
  mag_mapping_path: /path/to/mag_mapping.tsv
```

---

### output

Output directory configuration.

| Parameter | Description |
|-----------|-------------|
| `root_dir` | Root directory for all output files. Subdirectories will be created for each analysis step. |

**Example:**

```yaml
output:
  root_dir: ./alleleflux_output
```

---

### analysis

Core analysis settings.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `data_type` | `longitudinal` | Type of analysis: `single` (one timepoint) or `longitudinal` (multiple timepoints). |
| `allele_analysis_only` | `false` | If true, only run allele frequency analysis without statistical tests. |
| `use_lmm` | `true` | Enable Linear Mixed Models (LMM) for repeated measures/longitudinal data. Best for accounting for subject-level variation. |
| `use_significance_tests` | `true` | Enable two-sample (t-test, Mann-Whitney) and single-sample statistical tests. Best for simple comparisons. |
| `use_cmh` | `true` | Enable Cochran-Mantel-Haenszel tests for stratified categorical analysis. Best for detecting consistent directional changes. |
| `timepoints_combinations` | Required | List of timepoint combinations to analyze (see below). |
| `groups_combinations` | Required | List of group pairs to compare (see below). |

:::{seealso}
For detailed information about statistical tests and score calculations, see [Statistical Tests Reference](statistical_tests.md).
:::

**Timepoints Configuration:**

For longitudinal analysis, specify pairs of timepoints and a **focus timepoint**:

```yaml
analysis:
  data_type: longitudinal
  timepoints_combinations:
    - timepoint: [pre, post]
      focus: post      # The later/derived timepoint
    - timepoint: [pre, mid]
      focus: mid
```

**Understanding the Focus Timepoint:**

The focus timepoint represents the **derived** or **later** state in evolutionary comparisons:

- **For dN/dS analysis**: The focus timepoint is treated as the "derived" (Time 2) state, while the other timepoint is "ancestral" (Time 1). AlleleFlux calculates evolutionary changes in the direction: ancestral → derived.
- **For CMH scores**: The score measures differential significance relative to the focus timepoint (sites significant at focus but not at the other timepoint).
- **Selection guideline**: Always choose the **later** or **endpoint** timepoint as focus.
- **Default behavior**: If not specified, defaults to the second timepoint in the list.

**Examples:**

```yaml
# Typical longitudinal study: Day 0 → Day 30
timepoints_combinations:
  - timepoint: [day0, day30]
    focus: day30        # day30 is derived, day0 is ancestral

# Treatment study: Baseline → Post-treatment
timepoints_combinations:
  - timepoint: [baseline, post_treatment]
    focus: post_treatment  # Post is derived state
```

For single timepoint analysis:

```yaml
analysis:
  data_type: single
  timepoints_combinations:
    - timepoint: [baseline]  # No focus needed for single timepoint
```

**Groups Configuration:**

Specify pairs of groups to compare:

```yaml
analysis:
  groups_combinations:
    - [treatment, control]
    - [high_fat, standard]
```

---

### quality_control

Parameters for filtering samples and positions.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_sample_num` | `4` | Minimum number of samples per group required for statistical tests. MAGs with fewer valid samples are marked as ineligible. |
| `breadth_threshold` | `0.1` | Minimum breadth of coverage (fraction of genome with ≥1x coverage). Range: 0.0-1.0. |
| `coverage_threshold` | `1.0` | Minimum average coverage depth required. Samples below this are excluded. |
| `disable_zero_diff_filtering` | `false` | If true, keep positions where allele frequencies do not change. By default, constant positions are filtered out. |

**Example:**

```yaml
quality_control:
  min_sample_num: 4
  breadth_threshold: 0.1
  coverage_threshold: 1.0
  disable_zero_diff_filtering: false
```

---

### profiling

Parameters for BAM file processing during profiling.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ignore_orphans` | `true` | Ignore orphan reads (unpaired reads). Set to `false` to include unpaired reads. |
| `min_base_quality` | `30` | Minimum Phred base quality score to include a base in the pileup. |
| `min_mapping_quality` | `2` | Minimum mapping quality (MAPQ) score to include a read. |
| `ignore_overlaps` | `true` | Ignore overlapping segments of read pairs to avoid double-counting. |

**Example:**

```yaml
profiling:
  ignore_orphans: true
  min_base_quality: 30
  min_mapping_quality: 2
  ignore_overlaps: true
```

:::{note}
Higher `min_base_quality` values (e.g., 30) reduce sequencing errors but may also reduce coverage. For high-quality data, the default of 30 is recommended.
:::

---

### statistics

Parameters for statistical testing.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `filter_type` | `t-test` | Type of initial filter for preprocessing positions. |
| `preprocess_between_groups` | `true` | Enable preprocessing for between-group comparisons. |
| `preprocess_within_groups` | `true` | Enable preprocessing for within-group comparisons. |
| `max_zero_count` | `4` | Maximum number of zero-frequency samples allowed per position in preprocessing. |
| `p_value_threshold` | `0.05` | Significance threshold (alpha) for statistical tests. |
| `fdr_group_by_mag_id` | `false` | If true, apply FDR correction within each MAG. If false, apply across all positions. |
| `min_positions_after_preprocess` | `1` | Minimum number of positions required after preprocessing to proceed with analysis. |

**Example:**

```yaml
statistics:
  filter_type: t-test
  preprocess_between_groups: true
  preprocess_within_groups: true
  max_zero_count: 4
  p_value_threshold: 0.05
  fdr_group_by_mag_id: false
  min_positions_after_preprocess: 1
```

---

### dnds

Parameters for dN/dS (synonymous/non-synonymous) ratio calculations.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `p_value_column` | `q_value` | Column name to use for significance in dN/dS calculations. |
| `dn_ds_test_type` | `two_sample_unpaired_tTest` | Type of statistical test to use for dN/dS analysis. |

**Example:**

```yaml
dnds:
  p_value_column: q_value
  dn_ds_test_type: two_sample_unpaired_tTest
```

---

### resources

Computational resource allocation for cluster execution.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `threads_per_job` | `16` | Number of CPU threads allocated to each job. |
| `mem_per_job` | `8G` | Memory allocation per job. Formats: `8G`, `16GB`, `8192M`. |
| `time` | `24:00:00` | Maximum wall time per job in HH:MM:SS format. |

**Example:**

```yaml
resources:
  threads_per_job: 16
  mem_per_job: 8G
  time: '24:00:00'
```

---

## Complete Configuration Example

```yaml
run_name: diet_microbiome_study

input:
  fasta_path: /data/mags/combined_mags.fasta
  prodigal_path: /data/mags/prodigal_genes.fna
  metadata_path: /data/metadata/samples_with_bam.tsv
  gtdb_path: /data/taxonomy/gtdbtk.bac120.summary.tsv
  mag_mapping_path: /data/mags/mag_mapping.tsv

output:
  root_dir: ./results

log_level: INFO

analysis:
  data_type: longitudinal
  allele_analysis_only: false
  use_lmm: true
  use_significance_tests: true
  use_cmh: true
  timepoints_combinations:
    - timepoint: [pre, post]
      focus: post
  groups_combinations:
    - [high_fat, control]

quality_control:
  min_sample_num: 4
  breadth_threshold: 0.1
  coverage_threshold: 1.0
  disable_zero_diff_filtering: false

profiling:
  ignore_orphans: true
  min_base_quality: 30
  min_mapping_quality: 2
  ignore_overlaps: true

statistics:
  filter_type: t-test
  preprocess_between_groups: true
  preprocess_within_groups: true
  max_zero_count: 4
  p_value_threshold: 0.05
  fdr_group_by_mag_id: false
  min_positions_after_preprocess: 1

dnds:
  p_value_column: q_value
  dn_ds_test_type: two_sample_unpaired_tTest

resources:
  threads_per_job: 16
  mem_per_job: 8G
  time: '24:00:00'
```

## Quick Tips

- **breadth_threshold**: Start with `0.1` (10% coverage); increase for high-coverage data
- **min_sample_num**: Minimum `4` samples per group for robust inference
- **min_base_quality**: Keep at `30` for Illumina; lower to `20` for older data
- **Resource allocation**: Adjust `threads_per_job` and `mem_per_job` based on MAG sizes

For worked examples, see [Use Cases](../examples/use_cases.md).
