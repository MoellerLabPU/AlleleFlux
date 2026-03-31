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
- `total_coverage`: Total read coverage
- `A`, `C`, `G`, `T`, `N`: Base counts
- `mapq_scores`: MAPQ scores for reads
- `gene_id`: Overlapping gene identifier (if any)

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
| `--coverage-threshold` | 1.0 | Minimum average coverage depth. |
| `--data-type` | `longitudinal` | Analysis type: `single` or `longitudinal`. |

#### Output

Creates `{mag_id}_QC.tsv` with QC results including `breadth_threshold_passed` column.

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
| `--metadata` | Enhanced metadata file with `sample_profile_dir` column. |
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
  --profile-dir DIR --metadata META --group GROUP \
  --timepoint TP --output DIR [options]
```

#### Required Arguments

| Argument | Description |
|----------|-------------|
| `--significant-sites` | Path to significant sites table (from p-value summary). |
| `--profile-dir` | Directory containing sample profile subdirectories. |
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
