# CLI Reference

This section provides a complete reference for all AlleleFlux command-line tools.

## Main CLI

The main `alleleflux` command provides access to workflow execution and configuration:

```bash
alleleflux --help
```

Subcommands:
: - `run`: Execute the AlleleFlux workflow
  - `init`: Create a new configuration file interactively
  - `info`: Display installation information
  - `tools`: List all available CLI tools

### alleleflux run

Execute the complete AlleleFlux pipeline.

```bash
alleleflux run --config CONFIG_FILE [OPTIONS]
```

**Required Options:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Option
     - Description
   * - ``-c, --config``
     - Path to the AlleleFlux configuration file (YAML)
```

**Optional Options:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``-w, --working-dir``
     - ``.``
     - Working directory for the workflow
   * - ``-j, --jobs``
     - None
     - Max concurrent jobs (local execution only)
   * - ``-t, --threads``
     - All CPUs
     - Total threads available (local execution only)
   * - ``-m, --memory``
     - None
     - Total memory available (e.g., '64G', '128GB')
   * - ``-p, --profile``
     - None
     - Snakemake profile directory for cluster execution
   * - ``-n, --dry-run``
     - False
     - Perform a dry run without executing jobs
   * - ``--unlock``
     - False
     - Unlock the working directory (removes stale locks)
   * - ``--snakemake-args``
     - None
     - Additional Snakemake arguments (quoted string)
```

**Examples:**

```bash
# Run with a config file
alleleflux run --config config.yml

# Run with limited resources
alleleflux run --config config.yml --threads 16 --memory 64G

# Dry run to see what would be executed
alleleflux run --config config.yml --dry-run

# Run with SLURM profile
alleleflux run --config config.yml --profile slurm_profile/

# Force rerun all jobs
alleleflux run --config config.yml -- --forceall --reason
```

### alleleflux init

Create a new configuration file interactively.

```bash
alleleflux init [OPTIONS]
```

**Options:**

```{eval-rst}
.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Option
     - Default
     - Description
   * - ``--template``
     - False
     - Print a template config to stdout
   * - ``--output``
     - ``alleleflux_config.yml``
     - Output config file path
```

\---

## Analysis Scripts

### alleleflux-profile

Profile MAGs using alignment files. Processes BAM files to extract base-level coverage and allele information.

```bash
alleleflux-profile --bam_path BAM --fasta_path FASTA --prodigal_fasta GENES --mag_mapping_file MAPPING --output_dir DIR [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--bam_path``
     - Path to sorted BAM file
   * - ``--fasta_path``
     - Path to reference FASTA file used for alignment
   * - ``--prodigal_fasta``
     - Path to Prodigal predicted genes (DNA FASTA)
   * - ``--mag_mapping_file``
     - Tab-separated file mapping contig names to MAG IDs (columns: ``contig_name``, ``mag_id``)
   * - ``--output_dir``
     - Path to output directory
```

**Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--cpus``
     - All CPUs
     - Number of processors to use
   * - ``--sampleID``
     - From BAM filename
     - Sample identifier (overrides automatic extraction)
   * - ``--no-ignore-orphans``
     - False
     - Include reads without a properly paired mate
   * - ``--min-base-quality``
     - 30
     - Minimum base quality score to include a base
   * - ``--min-mapping-quality``
     - 2
     - Minimum mapping quality score to include a read
   * - ``--no-ignore-overlaps``
     - False
     - Do not ignore overlapping read segments (may double-count)
   * - ``--log-level``
     - INFO
     - Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
```

**Output:**

Creates `{output_dir}/{sampleID}/{sampleID}_{mag_name}_profiled.tsv.gz` files with columns:

- `contig`: Contig identifier
- `position`: 0-based genomic position
- `ref_base`: Reference base at this position
- `total_coverage`: Total read coverage
- `A`, `C`, `G`, `T`, `N`: Base counts
- `mapq_scores`: MAPQ scores for reads at this position
- `gene_id`: Overlapping gene identifier (if any)

### alleleflux-allele-freq

Analyze allele frequencies across samples.

```bash
alleleflux-allele-freq --mag_id MAG_ID --mag_metadata_file METADATA --fasta FASTA --output_dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Description
   * - ``--mag_id``
     - MAG identifier to process
   * - ``--mag_metadata_file``
     - Path to MAG metadata file
   * - ``--fasta``
     - Path to reference FASTA file
   * - ``--breadth_threshold``
     - Minimum breadth of coverage (default: 0.1)
   * - ``--data_type``
     - Analysis type: ``single`` or ``longitudinal``
   * - ``--output_dir``
     - Path to output directory
   * - ``--disable_zero_diff_filtering``
     - Disable filtering of constant positions
```

### alleleflux-scores

Calculate parallelism and divergence scores.

```bash
alleleflux-scores --significance_df INPUT --test_type TEST --mag_id MAG --output_dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Argument
     - Description
   * - ``--significance_df``
     - Path to significance test results
   * - ``--test_type``
     - Type of test: ``two_sample_paired``, ``two_sample_unpaired``, ``single_sample``
   * - ``--mag_id``
     - MAG identifier
   * - ``--p_value_threshold``
     - Significance threshold (default: 0.05)
   * - ``--output_dir``
     - Path to output directory
```

### alleleflux-cmh-scores

Calculate CMH test scores for a MAG.

```bash
alleleflux-cmh-scores --cmh_df CMH_RESULTS --mag_id MAG --output_dir DIR [OPTIONS]
```

\---

## Preprocessing Scripts

### alleleflux-metadata

Generate MAG metadata files from sample profiles.

```bash
alleleflux-metadata --metadata_file INPUT --profiles_dir DIR --mag_id MAG --output_dir DIR [OPTIONS]
```

### alleleflux-qc

Perform quality control on MAG profiles.

```bash
alleleflux-qc --root_dir PROFILES --fasta FASTA --mag_mapping_file MAPPING --mag_id MAG --output_dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--breadth_threshold``
     - 0.1
     - Minimum breadth of coverage (0-1)
   * - ``--coverage_threshold``
     - 1.0
     - Minimum average coverage depth
   * - ``--data_type``
     - longitudinal
     - Analysis type: ``single`` or ``longitudinal``
```

### alleleflux-eligibility

Generate eligibility tables for MAGs.

```bash
alleleflux-eligibility --qc_dir QC_DIR --min_sample_num N --output_file OUTPUT [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--min_sample_num``
     - 4
     - Minimum number of samples required for statistical tests
   * - ``--data_type``
     - longitudinal
     - Analysis type: ``single`` or ``longitudinal``
```

### alleleflux-preprocessing-eligibility

Check preprocessing eligibility for MAGs.

```bash
alleleflux-preprocessing-eligibility --input_dir DIR --output_file OUTPUT [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--input_dir``
     - Required
     - Directory containing preprocessed data
   * - ``--output_file``
     - Required
     - Output eligibility file path
   * - ``--min_positions_after_preprocess``
     - 1
     - Minimum positions required after preprocessing
```

### alleleflux-p-value-summary

Summarize p-values across statistical tests.

```bash
alleleflux-p-value-summary --input_dir DIR --output_file OUTPUT [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--input_dir``
     - Required
     - Directory containing statistical test results
   * - ``--output_file``
     - Required
     - Output summary file path
   * - ``--fdr_group_by_mag_id``
     - False
     - Apply FDR correction per MAG
   * - ``--p_value_threshold``
     - 0.05
     - P-value threshold for filtering
```

\---

## Statistics Scripts

### alleleflux-cmh

Run Cochran-Mantel-Haenszel tests for stratified analysis.

```bash
alleleflux-cmh --input_df INPUT --mag_id MAG --output_dir DIR [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--input_df``
     - Path to input allele frequency dataframe
   * - ``--mag_id``
     - MAG ID to process
   * - ``--output_dir``
     - Path to output directory
```

**Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--preprocessed_df``
     - None
     - Path to filtered dataframe for position filtering
   * - ``--min_sample_num``
     - 4
     - Minimum number of strata (replicates) required
   * - ``--data_type``
     - longitudinal
     - Analysis mode: ``single``, ``longitudinal``, or ``across_time``
   * - ``--group``
     - None
     - Group name for ``across_time`` mode
   * - ``--cpus``
     - All CPUs
     - Number of processors to use
```

**Output:**

Creates `{mag_id}_cmh.tsv.gz` with columns:

- `mag_id`: MAG identifier
- `contig`: Contig identifier
- `gene_id`: Gene identifier
- `position`: 0-based position
- `num_pairs`: Number of replicate pairs tested
- `p_value_CMH`: CMH test p-value
- `time`: Timepoint (for longitudinal data)
- `notes`: Error messages or warnings

### alleleflux-lmm

Run Linear Mixed Models for longitudinal analysis.

```bash
alleleflux-lmm --input_df INPUT --preprocessed_df PREPROCESSED --group GROUP --mag_id MAG --output_dir DIR [OPTIONS]
```

### alleleflux-single-sample

Perform single-sample statistical tests.

```bash
alleleflux-single-sample --input_df INPUT --group GROUP --mag_id MAG --output_dir DIR [OPTIONS]
```

### alleleflux-two-sample-paired

Perform paired two-sample tests.

```bash
alleleflux-two-sample-paired --input_df INPUT --mag_id MAG --output_dir DIR [OPTIONS]
```

### alleleflux-two-sample-unpaired

Perform unpaired two-sample tests.

```bash
alleleflux-two-sample-unpaired --input_df INPUT --mag_id MAG --output_dir DIR [OPTIONS]
```

\---

## Accessory Scripts

### alleleflux-create-mag-mapping

Create MAG mapping files from individual FASTA files.

```bash
alleleflux-create-mag-mapping --dir MAG_DIR --extension EXT --output-fasta COMBINED --output-mapping MAPPING [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--dir``
     - Directory containing individual MAG FASTA files
   * - ``--extension``
     - File extension of MAG files (e.g., ``fa``, ``fasta``)
   * - ``--output-fasta``
     - Path for combined output FASTA
   * - ``--output-mapping``
     - Path for contig-to-MAG mapping file
```

### alleleflux-add-bam-path

Add BAM file paths to metadata.

```bash
alleleflux-add-bam-path --metadata INPUT --output OUTPUT --bam-dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--metadata``
     - Required
     - Path to existing metadata file
   * - ``--output``
     - Required
     - Path to save updated metadata file
   * - ``--bam-dir``
     - ``.``
     - Directory containing BAM files
   * - ``--bam-extension``
     - ``.bam``
     - Extension of BAM files
   * - ``--drop-missing``
     - False
     - Drop samples without matching BAM files
```

### alleleflux-coverage-allele-stats

Calculate coverage and allele statistics for MAGs.

```bash
alleleflux-coverage-allele-stats --input_dir DIR --output_file OUTPUT [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--input_dir``
     - Required
     - Directory containing profile files
   * - ``--output_file``
     - Required
     - Output statistics file path
   * - ``--cpus``
     - All CPUs
     - Number of processors to use
```

**Output:**

Coverage and allele statistics per MAG including mean coverage, breadth, and allele diversity metrics.

### alleleflux-list-mags

List all MAG IDs found in a directory.

```bash
alleleflux-list-mags --input_dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--input_dir``
     - Required
     - Directory containing MAG files or profiles
   * - ``--output_file``
     - None
     - Optional output file (prints to stdout if not specified)
   * - ``--pattern``
     - ``*``
     - Glob pattern for file matching
```

### alleleflux-positions-qc

Perform quality control on genomic positions.

```bash
alleleflux-positions-qc --input_file FILE --output_file OUTPUT [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--input_file``
     - Required
     - Input position data file
   * - ``--output_file``
     - Required
     - Output QC results file
   * - ``--min_coverage``
     - 1
     - Minimum coverage threshold
   * - ``--min_samples``
     - 1
     - Minimum samples with data at position
```

### alleleflux-copy-profiles

Copy MAG profile files to a new location.

```bash
alleleflux-copy-profiles --source_dir DIR --dest_dir DIR [OPTIONS]
```

**Key Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--source_dir``
     - Required
     - Source directory containing profiles
   * - ``--dest_dir``
     - Required
     - Destination directory
   * - ``--mag_list``
     - None
     - Optional file with MAG IDs to copy (copies all if not specified)
   * - ``--symlink``
     - False
     - Create symlinks instead of copying
```

\---

## Visualization Scripts

AlleleFlux provides a complete visualization workflow. See {doc}`../usage/visualization_guide` for the complete workflow tutorial.

### alleleflux-prepare-metadata

Prepare and combine metadata tables for the visualization workflow.

```bash
alleleflux-prepare-metadata --metadata-in INPUT --metadata-out OUTPUT --base-profile-dir DIR [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--metadata-in``
     - Input metadata table (TSV)
   * - ``--metadata-out``
     - Output standardized metadata file (appends if exists)
   * - ``--base-profile-dir``
     - Base directory containing sample profile subdirectories
```

**Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--sample-col``
     - ``sample_id``
     - Column name for sample IDs in input
   * - ``--group-col``
     - ``group``
     - Column name for experimental groups
   * - ``--time-col``
     - ``time``
     - Column name for timepoints
   * - ``--day-col``
     - ``day``
     - Column name for day/order (optional)
   * - ``--replicate-col``
     - ``replicate``
     - Column name for replicates (optional)
   * - ``--subject-col``
     - ``subjectID``
     - Column name for subject IDs
```

**Output:**

Creates a standardized metadata file with columns: `sample_id`, `group`, `time`, `subjectID`, `sample_profile_dir`.

### alleleflux-terminal-nucleotide

Identify terminal nucleotides at significant genomic sites.

```bash
alleleflux-terminal-nucleotide --significant_sites SITES --profile_dir DIR --metadata META --group GROUP --timepoint TP --output DIR [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--significant_sites``
     - Path to significant sites table (from p-value summary)
   * - ``--profile_dir``
     - Directory containing sample profile subdirectories
   * - ``--metadata``
     - Sample metadata file
   * - ``--group``
     - Target group name for terminal nucleotide calculation
   * - ``--timepoint``
     - Target timepoint (typically endpoint)
   * - ``--output``
     - Output directory
```

**Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--p_value_column``
     - ``q_value``
     - Column for significance filtering: ``min_p_value`` or ``q_value``
   * - ``--p_value_threshold``
     - 0.05
     - Maximum p-value to include site
   * - ``--test_type``
     - ``two_sample_paired_tTest``
     - Test type to filter sites
   * - ``--group_filter``
     - None
     - Optional additional group filter
   * - ``--cpus``
     - All CPUs
     - Number of processors to use
   * - ``--log-level``
     - INFO
     - Logging level
```

**Output:**

- `{mag_id}/{mag_id}_terminal_nucleotides.tsv`: Terminal alleles per site
- `{mag_id}/{mag_id}_frequencies.tsv`: Full frequency data
- `terminal_nucleotide_analysis_summary.tsv`: Summary across MAGs

### alleleflux-track-alleles

Track anchor allele frequencies across all samples and timepoints.

```bash
alleleflux-track-alleles --mag-id MAG --anchor-file FILE --metadata META --output-dir DIR [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--mag-id``
     - MAG identifier to process
   * - ``--anchor-file``
     - Path to terminal nucleotides file (from terminal analysis)
   * - ``--metadata``
     - Enhanced metadata file with ``sample_profile_dir`` column
   * - ``--output-dir``
     - Output directory
```

**Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 30 45
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--anchor-column``
     - ``terminal_nucleotide_mean_freq``
     - Anchor column: ``terminal_nucleotide_mean_freq`` or ``terminal_nucleotide_majority_vote``
   * - ``--min-cov-per-site``
     - 0
     - Minimum coverage required per site
   * - ``--cpus``
     - All CPUs
     - Number of processors to use
```

**Output:**

- `{mag_id}_frequency_table.wide.tsv`: Sites × samples matrix
- `{mag_id}_frequency_table.long.tsv`: Tidy format for R/plotting

### alleleflux-plot-trajectories

Generate allele frequency trajectory visualizations.

```bash
alleleflux-plot-trajectories --input_file FILE [OPTIONS]
```

**Required Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Argument
     - Description
   * - ``--input_file``
     - Long-format frequency table (from track-alleles)
```

**Key Optional Arguments:**

```{eval-rst}
.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Argument
     - Default
     - Description
   * - ``--value_col``
     - ``min_p_value``
     - Column for ranking sites: ``min_p_value``, ``q_value``
   * - ``--n_sites_line``
     - 10
     - Number of top sites for line plots (or "all")
   * - ``--n_sites_dist``
     - all
     - Number of sites for box/violin plots
   * - ``--x_col``
     - ``time``
     - X-axis column: ``time`` or ``day``
   * - ``--x_order``
     - None
     - Custom x-axis order (space-separated values)
   * - ``--plot_types``
     - ``line``
     - Plot types: ``line``, ``box``, ``violin`` (space-separated)
   * - ``--per_site``
     - False
     - Generate individual plots per site
   * - ``--n_sites_per_site``
     - None
     - Number of sites for per-site plots
   * - ``--output_dir``
     - ``./plots``
     - Output directory
   * - ``--output_format``
     - ``png``
     - Format: ``png``, ``pdf``, ``svg``
   * - ``--group_by_replicate``
     - False
     - Aggregate trajectories by replicate
   * - ``--bin_width``
     - None
     - Day binning width (requires ``day`` column)
   * - ``--min_samples_per_bin``
     - 1
     - Minimum samples per time bin
   * - ``--line_alpha``
     - 0.8
     - Line transparency (0-1)
```

**Output:**

- `{mag_id}_line_plot.{format}`: Combined line trajectories
- `{mag_id}_box_plot.{format}`: Box plots by timepoint
- `{mag_id}_violin_plot.{format}`: Violin plots by timepoint
- `per_site/{contig}_{position}_{gene}_line.{format}`: Per-site plots (if enabled)

\---

## Evolution Scripts

### alleleflux-dnds-from-timepoints

Calculate dN/dS ratios from timepoint comparisons.

```bash
alleleflux-dnds-from-timepoints --input INPUT --output OUTPUT [OPTIONS]
```

\---

## Getting Help

For help with any command, use the `-h` or `--help` flag:

```bash
alleleflux --help
alleleflux run --help
alleleflux-profile --help
alleleflux-cmh --help
```
