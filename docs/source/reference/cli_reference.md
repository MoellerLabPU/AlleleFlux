# CLI Reference

AlleleFlux ships a main `alleleflux` entrypoint plus console scripts that power the Snakemake workflow. Most users only need `alleleflux run`; the other tools are available for advanced or ad‑hoc use.

## Main commands

### `alleleflux run` — execute the workflow

Execute the complete AlleleFlux pipeline with flexible resource control and scheduling options.

```bash
alleleflux run --config config.yml [options] [-- <extra snakemake args>]
```

#### Arguments

| Option | Default | Description |
|--------|---------|-------------|
| `-c, --config` | (required) | Path to the AlleleFlux configuration YAML. |
| `-w, --working-dir` | `.` | Working directory for Snakemake execution. |
| `-j, --jobs` | None | Max concurrent jobs (local only; ignored when `--profile` is set). |
| `-t, --threads` | None | Total threads available for local runs. |
| `-m, --memory` | None | Total memory for local runs (e.g., `64G`, `128GB`). |
| `-p, --profile` | None | Snakemake profile directory for cluster/HPC execution (e.g., `slurm_profile/`). |
| `-n, --dry-run` | False | Plan the DAG without running jobs. |
| `--unlock` | False | Unlock a previously crashed working directory. |
| `--snakemake-args` | None | Quoted string of extra Snakemake flags (alternative to `--`). |

#### Examples

```bash
# Run with a config file
alleleflux run --config config.yml

# Run with limited resources
alleleflux run --config config.yml --threads 16 --memory 64G

# Dry run to see what would be executed
alleleflux run --config config.yml --dry-run

# Run with SLURM profile
alleleflux run --config config.yml --profile slurm_profile/

# Force rerun all jobs with reasoning
alleleflux run --config config.yml -- --forceall --reason

# Run with specific working directory
alleleflux run --config config.yml --working-dir /path/to/workdir
```

#### Notes

- Pass additional Snakemake flags either after `--` or via `--snakemake-args`.
- When using `--profile`, job/thread/memory parameters are overridden by profile settings.
- See [Running the Workflow](../usage/running_workflow.md) for detailed scheduling instructions.

### `alleleflux init` — create a config

Create a new AlleleFlux configuration file interactively or from a template.

```bash
alleleflux init [--template] [--output alleleflux_config.yml]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--template` | False | Print a template config to stdout instead of interactive mode. |
| `--output` | `alleleflux_config.yml` | Output configuration file path. |

#### Examples

```bash
# Interactive mode (prompts for settings)
alleleflux init

# Print template to stdout
alleleflux init --template

# Interactive mode with custom output file
alleleflux init --output my_alleleflux_config.yml

# Save template to file
alleleflux init --template > my_template.yml
```

### `alleleflux info` — show install paths

Display version, package location, and Snakefile paths. Useful for debugging installation issues.

```bash
alleleflux info
```

### `alleleflux tools` — list console scripts

List all available console scripts grouped by functional category.

```bash
alleleflux tools [--category {Analysis,Preprocessing,Statistics,Evolution,Accessory,Visualization}]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--category` | None | Filter by category (optional). Lists all if not specified. |

#### Examples

```bash
# List all tools
alleleflux tools

# List only Analysis tools
alleleflux tools --category Analysis

# List Preprocessing tools
alleleflux tools --category Preprocessing
```

## Console scripts by stage

These are invoked automatically by the workflow but can be run manually for testing or custom tasks. Run any script with `--help` for full arguments.

---

## Analysis tools

### `alleleflux-profile` — profile BAM files into per-MAG allele tables

Extract base-level coverage and allele information from aligned BAM files.

```bash
alleleflux-profile --bam-path BAM --fasta-path FASTA --prodigal-fasta GENES \
  --mag-mapping-file MAPPING --output-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--bam-path` | Path to sorted BAM file. |
| `--fasta-path` | Path to reference FASTA file (must match BAM alignment reference). |
| `--prodigal-fasta` | Path to Prodigal predicted genes (DNA FASTA format). |
| `--mag-mapping-file` | Tab-separated file mapping contigs to MAG IDs (columns: `contig_name`, `mag_id`). |
| `--output-dir` | Output directory for profiles. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--cpus` | All available CPUs | Number of processors to use. |
| `--sampleID` | From BAM filename | Sample identifier (auto-extracted if not provided). |
| `--min-base-quality` | 30 | Minimum base quality score to include a base. |
| `--min-mapping-quality` | 2 | Minimum mapping quality score to include a read. |
| `--no-ignore-orphans` | False | Include reads without properly paired mate. |
| `--no-ignore-overlaps` | False | Do not ignore overlapping read segments (may double-count). |
| `--log-level` | INFO | Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). |

#### Output

Creates `{output_dir}/{sampleID}/{sampleID}_{mag_id}_profiled.tsv.gz` with columns:
- `contig`: Contig identifier
- `position`: 0-based genomic position
- `ref_base`: Reference base at position
- `total_coverage`: Read coverage, equal to `A + C + G + T`
- `A`, `C`, `G`, `T`: Base counts (ambiguous `N` calls are not counted; positions with no A/C/G/T read are not written)
- `gene_id`: Overlapping gene identifier (if any)

Profiles written by earlier versions carried an extra `N` column summed into `total_coverage`; they are not supported by this version. Re-profile.

#### Examples

```bash
# Basic profiling
alleleflux-profile --bam-path sample1.bam --fasta-path reference.fa \
  --prodigal-fasta genes.fna --mag-mapping-file mag_mapping.tsv \
  --output-dir profiles/

# With custom sample ID and resource limits
alleleflux-profile --bam-path sample1.bam --fasta-path reference.fa \
  --prodigal-fasta genes.fna --mag-mapping-file mag_mapping.tsv \
  --output-dir profiles/ --sampleID my_sample --cpus 8

# Stricter quality filtering
alleleflux-profile --bam-path sample1.bam --fasta-path reference.fa \
  --prodigal-fasta genes.fna --mag-mapping-file mag_mapping.tsv \
  --output-dir profiles/ --min-base-quality 35 --min-mapping-quality 10
```

---

### `alleleflux-allele-freq` — compute allele frequencies per MAG

Analyze allele frequencies across samples and timepoints.

```bash
alleleflux-allele-freq --mag-id MAG --mag-metadata-file METADATA \
  --fasta FASTA --output-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--mag-id` | MAG identifier to process. |
| `--mag-metadata-file` | Path to MAG metadata file (TSV with sample_id, file_path, group, time). |
| `--fasta` | Path to reference FASTA file. |
| `--output-dir` | Output directory for results. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--breadth-threshold` | 0.1 | Minimum breadth of coverage (0-1). |
| `--data-type` | `longitudinal` | Analysis type: `single` or `longitudinal`. |
| `--disable-zero-diff-filtering` | False | Keep constant positions (all samples same allele). |

#### Output

Creates `{mag_id}_allele_freq.tsv.gz` with allele frequency data per position/sample.

---

### `alleleflux-scores` — calculate parallelism and divergence scores

Derive MAG-level scores from statistical test results.

```bash
alleleflux-scores --rootDir DIR --output-dir DIR [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--rootDir` | Required | Directory containing `*_metadata.tsv` files. |
| `--output-dir` | Required | Output directory. |
| `--cpus` | All available | Number of processors. |

#### Examples

```bash
# Score all MAGs
alleleflux-scores --rootDir metadata/ --output-dir scores/

# With custom CPU count
alleleflux-scores --rootDir metadata/ --output-dir scores/ --cpus 16
```

---

### `alleleflux-cmh-scores` — CMH-specific score aggregation

Calculate CMH test scores for a MAG. See also: `alleleflux-cmh` for running CMH tests.

```bash
alleleflux-cmh-scores --cmh-df INPUT --mag-id MAG --output-dir DIR [options]
```

---

## Preprocessing tools

### `alleleflux-metadata` — build MAG metadata from profiles

Generate MAG metadata files from sample profiles and sample sheet.

```bash
alleleflux-metadata --metadata-file INPUT --profiles-dir DIR \
  --mag-id MAG --output-dir DIR [options]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `--metadata-file` | Input sample metadata file (CSV/TSV). |
| `--profiles-dir` | Directory containing profile files. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

---

### `alleleflux-qc` — quality control on profiles

Perform coverage and breadth QC on MAG profiles.

```bash
alleleflux-qc --root-dir PROFILES --mag-id MAG --output-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--root-dir` | Directory containing profile files. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--fasta` | None | Path to reference FASTA (optional). |
| `--mag-mapping-file` | None | Contig-to-MAG mapping file (optional). |
| `--breadth-threshold` | 0.1 | Minimum breadth of coverage (0-1). |
| `--coverage-threshold` | 1.0 | Minimum `average_coverage_including_zeros` (summed depth / genome size). |
| `--data-type` | `longitudinal` | Analysis type: `single` or `longitudinal`. |

#### Output

Creates `{mag_id}_QC.tsv`, one row per sample, with the sample metadata, `genome_size`, and:

- `breadth`, `breadth_threshold_passed`, `breadth_fail_reason`: fraction of the genome with at least one read, and the breadth check
- `average_coverage`: mean depth over the **covered** positions only (summed depth / number of positions with at least one read). This is the reference the per-site depth filter divides by (see `quality_control.site_depth_filter` in the config); it is computed once here and read from this table downstream.
- `average_coverage_including_zeros`: summed depth / genome size, with uncovered positions counted as zero. This is the number `--coverage-threshold` tests (`coverage_threshold_passed`, `coverage_fail_reason`).
- `median_coverage`, `median_coverage_including_zeros`, `coverage_std`, `coverage_std_including_zeros`, `length_weighted_coverage`: further depth statistics (the `_including_zeros` ones treat uncovered positions as zero)

The two averages answer different questions: a genome with a handful of very deep positions can have a high `average_coverage` and a near-zero `average_coverage_including_zeros`, and it is the latter that decides whether the sample is kept.

#### Examples

```bash
# Basic QC
alleleflux-qc --root-dir profiles/ --mag-id MAG000001 --output-dir qc/

# Custom thresholds
alleleflux-qc --root-dir profiles/ --mag-id MAG000001 --output-dir qc/ \
  --breadth-threshold 0.2 --coverage-threshold 5.0
```

---

### `alleleflux-eligibility` — generate MAG eligibility tables

Create eligibility tables for statistical tests based on QC results.

```bash
alleleflux-eligibility --qc-dir QC_DIR --output-file OUTPUT [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--qc-dir` | Required | Directory containing QC files. |
| `--output-file` | Required | Output eligibility file path. |
| `--min-sample-num` | 4 | Minimum number of samples required. |
| `--data-type` | `longitudinal` | Analysis type: `single` or `longitudinal`. |

#### Output

Creates eligibility table with columns:
- `mag_id`: MAG identifier
- `unpaired_test_eligible`: Eligible for unpaired tests
- `paired_test_eligible`: Eligible for paired tests
- `single_sample_eligible_*`: Per-group single-sample eligibility

---

## Statistical test tools

### `alleleflux-cmh` — Cochran-Mantel-Haenszel stratified test

Run CMH tests for stratified allele frequency analysis (typically stratified by replicate).

```bash
alleleflux-cmh --input-df INPUT --mag-id MAG --output-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--input-df` | Path to input allele frequency dataframe. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--preprocessed-df` | None | Path to filtered dataframe for position filtering. |
| `--min-sample-num` | 4 | Minimum number of strata (replicates) required. |
| `--data-type` | `longitudinal` | Analysis mode: `single`, `longitudinal`, or `across_time`. |
| `--group` | None | Group name for `across_time` mode. |
| `--cpus` | All available | Number of processors. |

#### Output

Creates `{mag_id}_cmh.tsv.gz` with columns:
- `mag_id`: MAG identifier
- `contig`: Contig identifier
- `gene_id`: Gene identifier
- `position`: 0-based position
- `num_pairs`: Number of replicate pairs tested
- `p_value_CMH`: CMH test p-value
- `time`: Timepoint (for longitudinal data)
- `notes`: Error messages or warnings

#### Examples

```bash
# Basic CMH test
alleleflux-cmh --input-df allele_freq.tsv --mag-id MAG000001 --output-dir cmh_results/

# Across timepoints mode
alleleflux-cmh --input-df allele_freq.tsv --mag-id MAG000001 \
  --output-dir cmh_results/ --data-type across_time --group fat

# With preprocessing filter
alleleflux-cmh --input-df allele_freq.tsv --mag-id MAG000001 \
  --output-dir cmh_results/ --preprocessed-df preproc.tsv --cpus 16
```

---

### `alleleflux-lmm` — Linear mixed models for longitudinal analysis

Run LMM tests for longitudinal data with mixed effects.

```bash
alleleflux-lmm --input-df INPUT --preprocessed-df PREPROCESSED \
  --group GROUP --mag-id MAG --output-dir DIR [options]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `--input-df` | Path to input allele frequency dataframe. |
| `--preprocessed-df` | Path to filtered dataframe. |
| `--group` | Group name to analyze. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

---

### `alleleflux-two-sample-unpaired` — unpaired two-sample tests

Perform unpaired Mann-Whitney U tests comparing two groups.

```bash
alleleflux-two-sample-unpaired --input-df INPUT --mag-id MAG \
  --output-dir DIR [options]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `--input-df` | Path to input allele frequency dataframe. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

---

### `alleleflux-two-sample-paired` — paired two-sample tests

Perform paired Wilcoxon signed-rank tests on matched samples.

```bash
alleleflux-two-sample-paired --input-df INPUT --mag-id MAG \
  --output-dir DIR [options]
```

#### Arguments

| Argument | Description |
|----------|-------------|
| `--input-df` | Path to input allele frequency dataframe. |
| `--mag-id` | MAG ID to process. |
| `--output-dir` | Output directory. |

---

## Evolution tools

### `alleleflux-dnds-from-timepoints` — calculate dN/dS ratios

Compute dN/dS ratios from significant evolutionary sites.

```bash
alleleflux-dnds-from-timepoints --input INPUT --output OUTPUT [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--input` | Required | Path to input significant sites table. |
| `--output` | Required | Output dN/dS results file. |

See [dN/dS Analysis Guide](../usage/dnds_analysis.md) for detailed workflow.

---

## Accessory tools

### `alleleflux-create-mag-mapping` — generate MAG mapping and combined FASTA

Create contig-to-MAG mapping file and concatenate individual MAG FASTA files.

```bash
alleleflux-create-mag-mapping --dir MAG_DIR --extension EXT \
  --output-fasta COMBINED --output-mapping MAPPING [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--dir` | Directory containing individual MAG FASTA files. |
| `--extension` | File extension of MAG files (e.g., `fa`, `fasta`). |
| `--output-fasta` | Path for combined output FASTA. |
| `--output-mapping` | Path for contig-to-MAG mapping file (TSV). |

#### Output

- Combined FASTA: all contigs from all MAGs concatenated
- Mapping file: `contig_name\tmag_id` (tab-separated)

#### Examples

```bash
# Create mapping from directory of MAG FASTAs
alleleflux-create-mag-mapping --dir mags/ --extension fa \
  --output-fasta combined_reference.fa --output-mapping mag_mapping.tsv

# With different extension
alleleflux-create-mag-mapping --dir mags/ --extension fasta \
  --output-fasta reference.fasta --output-mapping mapping.tsv
```

---

### `alleleflux-add-bam-path` — add BAM file paths to metadata

Fill `bam_path` column in sample metadata by matching with BAM files.

```bash
alleleflux-add-bam-path --metadata INPUT --output OUTPUT \
  --bam-dir DIR [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--metadata` | Required | Path to input metadata file. |
| `--output` | Required | Path to save updated metadata. |
| `--bam-dir` | `.` | Directory containing BAM files. |
| `--bam-extension` | `.bam` | Extension of BAM files. |
| `--drop-missing` | False | Drop samples without matching BAM files. |

---

### `alleleflux-coverage-allele-stats` — compute coverage and allele statistics

Calculate coverage and allele statistics summary for all MAGs.

```bash
alleleflux-coverage-allele-stats --input-dir DIR --output-file OUTPUT [options]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--input-dir` | Required | Directory containing profile files. |
| `--output-file` | Required | Output statistics file path. |
| `--cpus` | All available | Number of processors. |

#### Output

Summary statistics per MAG: mean coverage, breadth, allele diversity metrics.

---

### `alleleflux-list-mags` — enumerate MAG IDs

List all unique MAG IDs from a directory of profile files.

```bash
alleleflux-list-mags --input-dir DIR [--output-file FILE]
```

#### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--input-dir` | Required | Directory containing MAG profile files. |
| `--output-file` | None | Optional output file (prints to stdout if not specified). |
| `--pattern` | `*` | Glob pattern for file matching. |

---

### Site-depth filter flags

Only `alleleflux-cache-allele-freq` reads raw profiles after QC, so it alone takes the per-site depth filter flags; the workflow passes them from the `quality_control.site_depth_filter` config block. Every other per-site consumer (pairwise ANI, baseline presence, dN/dS, the plotting tools) reads the cache that command writes, so the filter is applied exactly once.

| Argument | Default | Description |
|----------|---------|-------------|
| `--disable_site_depth_filter` | off | Do not apply the filter; profiles are used as-is. |
| `--site_depth` | 0 | Minimum depth at a position after counting; 0 = no floor (MIDAS2 uses 5). |
| `--site_ratio` | 3.0 | Drop a position when its depth divided by the sample's mean coverage over covered positions (QC's `average_coverage`) exceeds this. |
| `--min_allele_reads` | 0 | A base count must exceed this to count toward the position's depth; 0 = every read counts (MIDAS2 uses 2). |
| `--zero_low_counts` | off | MIDAS2 behaviour: at kept positions, base counts not exceeding `--min_allele_reads` are set to zero for everything downstream. Off = kept positions keep their original counts. |

### `alleleflux-site-filter-summary` — what the site-depth filter removed, run-wide

```bash
alleleflux-site-filter-summary --cache_dir RUN/allele_freq_cache --output site_depth_filter_summary.tsv
```

Run automatically by the workflow once every cache job has finished (local rule `site_depth_filter_summary`, outputs under `site_depth_filter/`); also usable by hand. Stacks every `*_site_depth_filter.tsv` sidecar the cache jobs wrote into one table (one row per sample × MAG × timepoint) and writes a per-MAG × timepoint rollup next to it as `<output stem>_by_mag.tsv` with `n_samples`, `mean_of_mean_coverage`, `mean_of_mean_coverage_after_filter`, `sites_removed_high_total`, `sites_covered_total` and `fraction_removed_high`. Fails loudly when the directory holds no sidecars.

---

### Additional accessory tools

- `alleleflux-positions-qc` — Position-level QC filtering
- `alleleflux-copy-profiles` — Copy or symlink profile files
- `alleleflux-single-sample` — Within-group single-sample test
- `alleleflux-preprocess-between-groups` — Position filtering between groups
- `alleleflux-preprocess-within-group` — Position filtering within groups
- `alleleflux-preprocessing-eligibility` — Aggregate preprocessing status
- `alleleflux-p-value-summary` — Summarize p-values across tests
- `alleleflux-outliers` — Flag outlier genes
- `alleleflux-taxa-scores` — Derive taxa-level scores
- `alleleflux-gene-scores` — Derive gene-level scores

---

## Visualization tools

### `alleleflux-plot-trajectories` — plot allele frequency trajectories

Generate allele frequency trajectory visualizations from tracked allele data.

```bash
alleleflux-plot-trajectories --input-file FILE [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--input-file` | Long-format frequency table from `alleleflux-track-alleles`. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--value-col` | `min_p_value` | Column for ranking sites: `min_p_value` or `q_value`. |
| `--n-sites-line` | 10 | Number of top sites for line plots (or `all`). |
| `--n-sites-dist` | `all` | Number of sites for box/violin plots. |
| `--x-col` | `time` | X-axis column: `time` or `day`. |
| `--x-order` | None | Custom x-axis order (space-separated values). |
| `--plot-types` | `line` | Plot types: `line`, `box`, `violin` (space-separated). |
| `--per-site` | False | Generate individual plots per site. |
| `--n-sites-per-site` | None | Number of sites for per-site plots. |
| `--output-dir` | `./plots` | Output directory. |
| `--output-format` | `png` | Format: `png`, `pdf`, `svg`. |
| `--group-by-replicate` | False | Aggregate trajectories by replicate. |
| `--bin-width` | None | Day binning width (requires `day` column). |
| `--min-samples-per-bin` | 1 | Minimum samples per time bin. |
| `--line-alpha` | 0.8 | Line transparency (0-1). |

#### Output

- `{mag_id}_line_plot.{format}`: Combined line trajectories
- `{mag_id}_box_plot.{format}`: Box plots by timepoint
- `{mag_id}_violin_plot.{format}`: Violin plots by timepoint
- `per_site/{contig}_{position}_{gene}_line.{format}`: Per-site plots (if enabled)

#### Examples

```bash
# Basic plotting
alleleflux-plot-trajectories --input-file tracked_alleles.tsv

# Multiple plot types with custom output
alleleflux-plot-trajectories --input-file tracked_alleles.tsv \
  --plot-types line box violin --output-dir results/plots/ \
  --output-format pdf

# Per-site plots for top 5 sites
alleleflux-plot-trajectories --input-file tracked_alleles.tsv \
  --per-site --n-sites-per-site 5 --output-format svg

# With binning and custom axis order
alleleflux-plot-trajectories --input-file tracked_alleles.tsv \
  --bin-width 7 --x-order "baseline week1 week2 week4 week8"
```

---

### `alleleflux-track-alleles` — track allele trajectories

Track anchor allele frequencies across all samples and timepoints.

```bash
alleleflux-track-alleles --mag-id MAG --anchor-file FILE \
  --metadata META --output-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--mag-id` | MAG identifier to process. |
| `--anchor-file` | Path to terminal nucleotides file (from `alleleflux-terminal-nucleotide`). |
| `--cache-dir` | The run's `allele_freq_cache/` directory (per-sample counts are read from the cache, never from raw profiles). |
| `--metadata` | Metadata file with `sample_id`, `group`, `time`, `subjectID`. |
| `--output-dir` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--anchor-column` | `terminal_nucleotide_mean_freq` | Anchor column to use for tracking. |
| `--min-cov-per-site` | 0 | Minimum coverage required per site. |
| `--cpus` | All available | Number of processors. |

#### Output

- `{mag_id}_frequency_table.wide.tsv`: Sites × samples matrix
- `{mag_id}_frequency_table.long.tsv`: Tidy format (for plotting)

---

### `alleleflux-prepare-metadata` — prepare metadata for visualization

Standardize and combine metadata tables for visualization workflows.

```bash
alleleflux-prepare-metadata --metadata-in INPUT --metadata-out OUTPUT \
  --base-profile-dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--metadata-in` | Input metadata table (TSV). |
| `--metadata-out` | Output standardized metadata file. |
| `--base-profile-dir` | Base directory containing sample profile subdirectories. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--sample-col` | `sample_id` | Column name for sample IDs. |
| `--group-col` | `group` | Column name for experimental groups. |
| `--time-col` | `time` | Column name for timepoints. |
| `--day-col` | `day` | Column name for day/order (optional). |
| `--replicate-col` | `replicate` | Column name for replicates (optional). |
| `--subject-col` | `subjectID` | Column name for subject IDs. |

#### Output

Standardized metadata with columns: `sample_id`, `group`, `time`, `subjectID`, `sample_profile_dir`.

---

### `alleleflux-terminal-nucleotide` — identify terminal nucleotides

Find dominant terminal alleles at significant genomic sites.

```bash
alleleflux-terminal-nucleotide --significant-sites SITES \
  --cache_dir DIR --metadata META --group GROUP \
  --timepoint TP --output DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--significant-sites` | Path to significant sites table (from p-value summary). |
| `--cache_dir` | The run's `allele_freq_cache/` directory (per-sample counts at the significant sites are read from the cache, never from raw profiles). |
| `--metadata` | Sample metadata file. |
| `--group` | Target group name for terminal nucleotide calculation. |
| `--timepoint` | Target timepoint (typically endpoint). |
| `--output` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--p-value-column` | `q_value` | Significance column: `min_p_value` or `q_value`. |
| `--p-value-threshold` | 0.05 | Maximum p-value to include site. |
| `--test-type` | `two_sample_paired_tTest` | Test type to filter sites. |
| `--group-filter` | None | Optional additional group filter. |
| `--cpus` | All available | Number of processors. |
| `--log-level` | INFO | Logging level. |

#### Output

- `{mag_id}/{mag_id}_terminal_nucleotides.tsv`: Terminal alleles per site
- `{mag_id}/{mag_id}_frequencies.tsv`: Full frequency data
- `terminal_nucleotide_analysis_summary.tsv`: Summary across MAGs

---

## Getting help

View detailed help for any tool:

```bash
# Main command help
alleleflux --help

# Subcommand help
alleleflux run --help
alleleflux init --help

# Console script help
alleleflux-profile --help
alleleflux-cmh --help
alleleflux-plot-trajectories --help
```

For configuration details, see [Configuration Reference](configuration.md). For how to run the workflow end to end, see [Running the Workflow](../usage/running_workflow.md).

## Strain turnover and baseline presence tools

See the [Strain Turnover and Baseline Presence guide](../usage/strain_turnover_analysis.md) for what these compute and how they chain.

### `alleleflux-pairwise-ani` — conANI / popANI between sample pairs of one MAG

```bash
alleleflux-pairwise-ani --mag MAG --cache_dir DIR --qc_files QC [QC ...] \
  --fasta FASTA --mag_mapping MAPPING --output_dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--mag` | MAG identifier. |
| `--cache_dir` | The run's `allele_freq_cache/` directory. Per-sample counts are read from the cache (QC-passing samples, site-depth-filtered), never from raw profiles. A QC-passing sample missing from an existing cache file stops the run; one whose timepoint has no cache file is skipped with a warning. |
| `--qc_files` | Per-MAG QC TSV(s), one per timepoint combination; samples that passed anywhere are compared. |
| `--fasta` | Reference FASTA. |
| `--mag_mapping` | MAG-to-contig mapping TSV (`mag_id`, `contig_id`). |
| `--output_dir` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--min_cov` | 5 | Reads required in **both** samples for a position to be compared; 1 disables the depth gate. |
| `--min_freq` | 0.05 | Minimum read fraction for an allele to count as present (popANI only). |
| `--fdr` | 1e-6 | Tolerated per-position probability that sequencing error explains an allele (popANI only). |
| `--min_base_quality` | 30 | Phred floor the profiles were built with; sets the error model. |
| `--pairs` | `within_subject` | `all`, `within_subject` (every same-mouse pair) or `transitions` (same-mouse pairs matching `--transitions`). |
| `--transitions` | — | `EARLIER:LATER` pairs for `--pairs transitions`, e.g. `pre:end`. |
| `--store_snp_locations` | `within_subject` | Which pairs also get per-position SNP rows: `none`, `within_subject`, `all`. |

#### Output

`{mag}_pairwise_ani.tsv` (one row per pair: `compared_bases_count`, `percent_genome_compared`, `consensus_SNPs`, `population_SNPs`, `conANI`, `popANI`, `min_cov`, and each sample's subject, group, time and replicate) and `{mag}_pairwise_ani_samples.tsv`.

---

### `alleleflux-strain-turnover` — per-mouse strain verdicts for one MAG

```bash
alleleflux-strain-turnover --mag MAG --pair_table TSV --output_dir DIR \
  --transitions EARLIER:LATER [EARLIER:LATER ...] [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--mag` | MAG identifier. |
| `--pair_table` | This MAG's `{mag}_pairwise_ani.tsv`. |
| `--output_dir` | Output directory. |
| `--transitions` | Transitions to call, e.g. `pre:end` or `5mo:10mo 5mo:22mo`. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--min_compared` | 0.1 | Minimum fraction of the genome compared for a verdict; below it both verdicts are `<NA>`. |
| `--pop_threshold` | 0.99999 | popANI below this = `strain_replacement`. |
| `--con_threshold` | 0.999 | conANI below this = `dominant_strain_change`. |

#### Output

`{mag}_strain_turnover.tsv` (one row per mouse × transition) and `{mag}_turnover_rollup.tsv` (counts per group × transition).

---

### `alleleflux-replacement-classification` — roll strain verdicts up to one row per MAG

```bash
alleleflux-replacement-classification --turnover_dir DIR --output_path TSV \
  [--min_voters 8] [--vote_rule majority] [--tie not_replaced] [--replicate_rule average]
```

| Argument | Description |
|----------|-------------|
| `--turnover_dir` | Directory of `{mag}_strain_turnover.tsv` files. |
| `--output_path` | Output TSV: one row per MAG × transition × metric, mouse block and replicate block side by side. |
| `--min_voters` | Fewest mice (or replicates) with a verdict needed for any yes/no call (default 8). Below it the flags are blank and `strain_status` reads `too_few_voters` or `no_voters`. |
| `--vote_rule` | When a checked key counts as replaced: `majority` (default), `any`, or `all`. |
| `--tie` | With `--vote_rule majority`, what exactly half the voters changing means for `strain_status`: `not_replaced` (default), `replaced`, or `unresolved`. |
| `--replicate_rule` | How the subjects of one replicate become its single vote in the replicate block: `average` (default: mean conANI / popANI of the subjects with a verdict, against the metric's threshold), `any`, or `majority`. |

---

### `alleleflux-rising-origin` — standing vs de novo for the alleles that rose

```bash
alleleflux-rising-origin --rising_dir DIR --output_dir DIR \
  [--classification TABLE [TABLE ...]] [--strain_metric METRIC] [--status_column COLUMN]
```

Reads every `*_baseline_presence_rising.tsv` in `--rising_dir` and counts, per comparison × group, the rising alleles in each origin category.

| Argument | Default | Description |
|----------|---------|-------------|
| `--rising_dir` | required | Folder of rising tables (`baseline_presence/`). |
| `--output_dir` | required | Output directory. |
| `--classification` | — | `replacement_classification.tsv` table(s), one per voter floor (`min_voters`). When given, de novo in a replaced MAG becomes `new_strain` and the summary is split by floor and strain class. |
| `--strain_metric` | `dominant_strain_change` | Classification metric (`dominant_strain_change` = conANI; `strain_replacement` = popANI, stricter). |
| `--status_column` | `strain_status_replicates` | Verdict column: replicates vote, or `strain_status` for mice. |
| `--cluster_window_bp` | `1000` | Cluster rule: two de novo alleles on the same contig that rose in the same group are neighbours when at most this many bp apart (any age); neighbours chain into blocks. |
| `--max_de_novo_neighbours` | `0` | Cluster rule: a de novo allele with more than this many de novo neighbours becomes a `block_event`, never an independent de novo mutation; 0 keeps only lone alleles as de novo. |
| `--turnover_dir` | — | Folder of `{mag}_strain_turnover.tsv` tables (`alleleflux-strain-turnover`): the per-mouse strain verdicts the carrier gate reads for each allele's carriers. |
| `--strain_gate` | `carrier` | Which strain verdict decides whether a de novo allele stands: `carrier` = the allele's own carriers (needs `--turnover_dir`); `group` = the across-cage vote (needs `--classification`); `both` = the stricter of the two (needs both inputs); `none` = no strain gate. See [Outputs](outputs.md#rising-origin). |
| `--min_determined_carriers` | `3` | Carrier gate: carriers with a per-mouse verdict needed before their strain is called; fewer → `strain_undetermined`. |

Outputs `rising_origin_alleles.tsv` and `rising_origin_summary.tsv`; see [Outputs](outputs.md).

---

### `alleleflux-baseline-presence` — was each significant allele present at baseline?

```bash
alleleflux-baseline-presence --run_dir DIR --comparison LABEL --test_type TYPE \
  --cache_dir DIR --metadata TSV --output_dir DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--run_dir` | Run root holding `p_value_summary/` and `significance_tests/`. |
| `--comparison` | Comparison label, e.g. `pre_end-fat_control`. |
| `--test_type` | Row filter within the summary, spelled exactly as the file does (`two_sample_paired_tTest`, `LMM_abs`, ...). |
| `--cache_dir` | The run's `allele_freq_cache/` directory. Verdicts are read from the cache (QC-passing samples only, site-filtered), never from raw profiles; a sample with no rows in a MAG's cache did not pass QC for that MAG and gets no rows. |
| `--metadata` | Sample metadata TSV (`sample_id`, `subjectID`, `group`, `time`[, `replicate`]). |
| `--output_dir` | Output directory. |

#### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--summary` | `two_sample_paired` | Summary family: `two_sample_paired`, `two_sample_unpaired`, `single_sample`, `lmm`, `lmm_across_time`. CMH is not offered (no per-base p-value). |
| `--threshold_column` | `q_value` | `q_value` or `min_p_value`. |
| `--threshold` | 0.05 | Significance cutoff. |
| `--baseline_scope` | `cohort` | What "absent at baseline" is judged against: `individual` = the carrier mice's own earlier samples; `own_group` = every earlier sample of the site's own group (within-group tests only: `single_sample`, `lmm_across_time`; an error with any other family); `comparison` = every earlier sample of both groups; `cohort` = every QC-passing sample of the run, any group, any timepoint, collected before the first carrier sample (needs `collection_date`, `YYYY-MM-DD`, for every sample in the metadata: a partly dated sheet is an error; with no dates at all "before" means an earlier timepoint label and a warning is logged). |
| `--turnover_dir` | — | Strain-turnover outputs; adds a `strain_background` column. |
| `--mags` | all | Restrict to these MAG ids. |
| `--min_cov` | 5 | Reads at a position before any verdict (1 = off). |
| `--min_freq` | 0.05 | Presence frequency floor. |
| `--fdr` | 1e-6 | Null-model false-discovery rate. |
| `--min_base_quality` | 30 | Base quality assumed by the null model. |
| `--disable_site_depth_filter`, `--site_depth`, `--site_ratio`, `--min_allele_reads`, `--zero_low_counts` | see below | The per-site depth filter, applied to each profile as it is loaded; a removed position reads as `not_covered` for that sample (see [Site-depth filter flags](#site-depth-filter-flags)). |
| `--cpus` | all | Worker processes (one profile load per MAG × sample). |

#### Output

`{comparison}_{family}_{statistic}_baseline_presence.tsv.gz` (long: one row per site × allele × sample) `..._summary.tsv` (one row per site × allele) and `..._rising.tsv` (one row per site × allele × group the allele rose in; see [Outputs](outputs.md)). Column names use the comparison's own timepoint labels.

