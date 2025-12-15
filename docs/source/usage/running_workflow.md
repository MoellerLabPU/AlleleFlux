# Running the Workflow

AlleleFlux uses Snakemake to manage the workflow. This guide explains how to run the workflow and understand its components.

## Workflow Overview

AlleleFlux is a unified Snakemake pipeline that performs:

1. **Profiling**: Process BAM files to extract allele frequencies for each MAG
2. **Quality Control**: Filter samples based on coverage breadth and depth
3. **Eligibility**: Determine which MAGs qualify for each statistical test
4. **Analysis**: Analyze allele frequencies across samples
5. **Statistical Testing**: Run appropriate tests based on your experimental design
6. **Scoring**: Calculate parallelism and divergence scores
7. **Outlier Detection**: Identify genes with exceptionally high scores

The pipeline automatically manages dependencies between these steps using Snakemake checkpoints.

## Running with the AlleleFlux CLI

The recommended way to run AlleleFlux is using the `alleleflux run` command:

```bash
alleleflux run --config config.yml
```

This command handles Snakemake invocation and resource management automatically.

**Common Options:**

```bash
# Specify threads and memory for local execution
alleleflux run --config config.yml --threads 16 --memory 64G

# Use a cluster profile for HPC execution
alleleflux run --config config.yml --profile slurm_profile/

# Dry run to see what would be executed
alleleflux run --config config.yml --dry-run

# Pass additional Snakemake arguments
alleleflux run --config config.yml -- --forceall --reason
```

## Running with Snakemake Directly

You can also run Snakemake directly for more control:

```bash
cd /path/to/AlleleFlux/alleleflux/smk_workflow/alleleflux_pipeline
snakemake --configfile /path/to/config.yml --profile profile/ --cores 16
```

## Available Statistical Tests

AlleleFlux supports several statistical testing approaches:

### Two-Sample Tests

These tests compare allele frequencies between groups:

- **Unpaired**: For comparing independent samples from different groups
- **Paired**: For comparing matched samples (e.g., before and after treatment)

### Single-Sample Tests

These tests analyze changes within a single group over time.

### Linear Mixed Models (LMM)

LMM tests account for complex experimental designs with fixed and random effects, particularly useful for longitudinal data.

### Cochran-Mantel-Haenszel (CMH) Tests

The CMH test is a stratified analysis of count data that:

- Tests for association between allele changes and conditions while controlling for confounding factors
- Provides position-by-position assessment of allele frequency changes
- Is especially powerful for detecting parallel evolutionary changes

To enable/disable specific tests, modify these settings in the `config.yml` file:

```yaml
analysis:
  use_significance_tests: true  # Enable/disable two-sample and single-sample tests
  use_lmm: true                # Enable/disable Linear Mixed Models
  use_cmh: true                # Enable/disable Cochran-Mantel-Haenszel tests
```

## Customizing the Workflow

You can customize the workflow by editing the `config.yml` file (see {doc}`input_preparation` for details). Key configuration options include:

```yaml
# Data type: "single" for a single timepoint or "longitudinal" for multiple timepoints
data_type: "longitudinal"

# Input files
input:
  bam_dir: "/path/to/bam_files"  # For backward compatibility
  fasta_path: "/path/to/reference.fa"
  prodigal_path: "/path/to/genes.fna"
  metadata_path: "/path/to/metadata.tsv"  # Must include bam_path column

# Output directory
output:
  root_dir: "/path/to/output"

# Quality control settings
quality_control:
  min_coverage_breadth: 0.5
  disable_zero_diff_filtering: false
  min_sample_num: 4
  breadth_threshold: 0.1

# Analysis settings
analysis:
  use_significance_tests: true
  use_lmm: true
  use_cmh: true
  significance_threshold: 0.05
```

## Advanced Usage

### Running on a Compute Cluster

AlleleFlux supports running on a compute cluster through Snakemake's cluster support. To run on a cluster:

1. Create a cluster profile for your system (see Snakemake documentation)
2. Run the workflow using your cluster profile:

```bash
# Using the AlleleFlux CLI with a cluster profile
alleleflux run --config config.yml --profile your_cluster_profile/

# Or using Snakemake directly
cd /path/to/AlleleFlux/alleleflux/smk_workflow/alleleflux_pipeline
snakemake --configfile /path/to/config.yml --profile your_cluster_profile/
```

Example cluster profiles are provided in the `profile/` directory. You can adapt these for your specific computing environment (e.g., SLURM, PBS, SGE).

### Output Files and Directories

The workflow generates several output directories:

```text
output/
├── profiles/               # Sample profiles
├── metadata/               # MAG metadata
├── eligibility/            # Eligibility tables
├── allele_analysis/        # Allele frequency analysis results
├── significance_tests/     # Statistical test results
│   ├── lmm/                # Linear Mixed Model results
│   ├── cmh/                # Cochran-Mantel-Haenszel test results
│   └── preprocessed_two_sample/  # Preprocessed data for two-sample tests
├── scores/                 # Parallelism and divergence scores
│   ├── per_MAG/            # Scores per MAG
│   └── processed/          # Processed scores (taxonomic and combined)
└── outliers/               # Outlier gene detection results
```

### Checkpoint Files

The workflow creates checkpoint files at various stages. If you need to restart a failed run, Snakemake will automatically pick up from the last successful checkpoint.

## Troubleshooting

If you encounter issues when running the workflow:

1. Check the Snakemake log files in the `logs/` directory
2. Ensure that all input files are in the correct format
3. Verify that you have sufficient resources (memory, CPU, disk space)
4. Check that all dependencies are installed correctly

### Common Issues

- **Error in rpy2 or R dependencies**: Ensure you have R installed and R packages required for CMH tests (e.g., stats)
- **Memory errors**: Increase the memory allocation in your Snakemake profile
- **Missing files**: Check paths in your config.yml file
