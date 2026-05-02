# Quickstart

Get started with AlleleFlux in minutes. This guide walks through configuration, execution, and output inspection.

## Prerequisites

Ensure you have:

- AlleleFlux installed (see [Installation](installation.md))
- Input files prepared:
  - BAM files (sorted and indexed, with `.bai` index files)
  - Reference FASTA (combined MAG contigs)
  - Prodigal gene predictions (`.fna` nucleotide output)
  - Sample metadata TSV (with `sample_id`, `bam_path`, `group`, `time` columns)
  - MAG mapping file (`contig_name` to `mag_id` assignments)

See [Input Preparation](../usage/input_preparation.md) for detailed file format specifications.

## 1. Initialize Configuration

The easiest way to create a configuration file is with the interactive wizard:

```bash
# Interactive configuration wizard (recommended)
alleleflux init
```

This prompts you for input file paths, analysis type, and parameters, then writes a ready-to-use `alleleflux_config.yml`.

Alternatively, generate a template to edit manually:

```bash
# Print the template to a file
alleleflux init --template > config.yml

# Or copy the template manually
cp $(python -c "import alleleflux; print(alleleflux.__path__[0])")/smk_workflow/config.template.yml config.yml
```

:::{tip}
`alleleflux init` is the recommended starting point for new users. It validates file paths and sets sensible defaults based on your analysis type (single vs. longitudinal).
:::

## 2. Edit Configuration

Open `config.yml` and set your file paths and analysis parameters. A minimal longitudinal configuration:

```yaml
run_name: "my_analysis"

input:
  fasta_path: "/data/combined_mags.fasta"
  prodigal_path: "/data/prodigal_genes.fna"
  metadata_path: "/data/sample_metadata.tsv"
  gtdb_path: "/data/gtdbtk.bac120.summary.tsv"
  mag_mapping_path: "/data/mag_mapping.tsv"

output:
  root_dir: "./alleleflux_output"

analysis:
  data_type: "longitudinal"
  use_significance_tests: true
  use_lmm: true
  use_cmh: true
  timepoints_combinations:
    - timepoint: ["pre", "post"]
      focus: "post"
  groups_combinations:
    - ["treatment", "control"]
```

For a single-timepoint study:

```yaml
analysis:
  data_type: "single"
  timepoints_combinations:
    - timepoint: ["baseline"]
  groups_combinations:
    - ["disease", "healthy"]
```

See [Configuration Reference](../reference/configuration.md) for all options.

## 3. Run the Pipeline

```bash
# Run locally with 16 threads
alleleflux run --config config.yml --threads 16

# Dry run to preview the execution plan
alleleflux run --config config.yml --dry-run

# Run with memory limit
alleleflux run --config config.yml --threads 16 --memory 64G
```

For HPC clusters with SLURM:

```bash
# Copy the bundled SLURM profile
cp -r $(python -c "import alleleflux; print(alleleflux.__path__[0])")/smk_workflow/slurm_profile ./

# Run with SLURM job submission
alleleflux run --config config.yml --profile ./slurm_profile
```

## 4. Examine Output

Results are organized under the output directory:

```text
alleleflux_output/
└── longitudinal/                # Subdirectory matches data_type
    ├── profiles/                # Per-sample allele frequency profiles
    ├── inputMetadata/           # Per-MAG metadata tables
    ├── QC/                      # Quality control metrics per MAG
    ├── eligibility_table_*.tsv  # Which MAGs qualify for each test
    ├── allele_analysis/         # Allele frequency analysis results
    ├── significance_tests/      # Statistical test results
    │   ├── two_sample_paired_*/ #   Paired t-test & Wilcoxon
    │   ├── two_sample_unpaired_*/ # Unpaired t-test & Mann-Whitney
    │   ├── single_sample_*/     #   Single-sample t-test
    │   ├── lmm_*/               #   Linear mixed models
    │   └── cmh_*/               #   CMH tests
    ├── scores/                  # Parallelism & divergence scores
    │   ├── intermediate/        #   Per-MAG raw scores
    │   └── processed/           #   Combined & taxonomic aggregations
    │       ├── combined/        #     MAG-level score tables
    │       └── gene_scores_*/   #     Gene-level score tables
    ├── outlier_genes/           # High-scoring genes (selection targets)
    └── dnds_analysis/           # dN/dS evolutionary rate analysis
```

**Key output files to check first:**

```bash
OUT=alleleflux_output/longitudinal

# 1. Check MAG eligibility
column -t -s $'\t' $OUT/eligibility_table_pre_post-treatment_control.tsv | head

# 2. View MAG-level scores (parallelism and divergence)
column -t -s $'\t' \
  $OUT/scores/processed/combined/MAG/scores_two_sample_paired-pre_post-treatment_control-MAGs.tsv

# 3. Inspect outlier genes for a specific MAG
zcat $OUT/outlier_genes/pre_post-treatment_control/*_outlier_genes.tsv.gz | head

# 4. View dN/dS results
zcat $OUT/dnds_analysis/pre_post-treatment_control/*_gene_summary_ng86.tsv.gz | head
```

See [Output Files Reference](../reference/outputs.md) for complete file format specifications.

## Understanding the Output

After a run completes, here is a recommended order for inspecting results:

1. **Eligibility tables** (`eligibility_table_*.tsv`): Start here. These tables show which MAGs passed quality control and are eligible for each statistical test. If a MAG you expect is missing from downstream results, check its eligibility status and QC metrics first.

2. **Scores** (`scores/processed/combined/`): The MAG-level score tables rank populations by parallelism and divergence. High parallelism scores indicate allele frequency changes occurring consistently across replicates -- a hallmark of selection. Sort by `parallelism_score` descending to find the most interesting MAGs.

3. **Outlier genes** (`outlier_genes/`): These are genes with scores significantly above the genome-wide background, representing candidate targets of selection. Each file contains gene IDs, scores, and functional annotations (if GTDB taxonomy was provided).

4. **Statistical tests** (`significance_tests/`): Per-position p-values from the various tests (paired/unpaired, LMM, CMH). These are primarily consumed by the scoring step but can be inspected directly for specific positions of interest.

5. **dN/dS analysis** (`dnds_analysis/`): Evolutionary rate estimates per gene, useful for distinguishing positive selection (dN/dS > 1) from neutral drift.

## 5. Run Individual Tools

Each analysis step is also available as a standalone command:

```bash
# Create MAG mapping from individual FASTA files
alleleflux-create-mag-mapping \
    --dir mag_fastas/ --extension fa \
    --output-fasta combined.fasta --output-mapping mapping.tsv

# Profile a single BAM file
alleleflux-profile \
    --bam_path sample.bam \
    --fasta_path reference.fa \
    --prodigal_fasta genes.fna \
    --mag_mapping_file mapping.tsv \
    --output_dir profiles/

# Quality control for a single MAG
alleleflux-qc \
    --root_dir profiles/ \
    --fasta reference.fa \
    --mag_mapping_file mapping.tsv \
    --mag_id Bacteroides_001 \
    --output_dir qc/ \
    --breadth_threshold 0.1

# Calculate dN/dS for a MAG
alleleflux-dnds-from-timepoints \
    --input p_value_summary.tsv \
    --output dnds/ \
    --mag_id Bacteroides_001 \
    --profiles_dir profiles/ \
    --prodigal_fasta genes.fna \
    --fasta reference.fa \
    --ancestral_timepoint pre \
    --derived_timepoint post
```

For complete CLI documentation: `alleleflux-<tool> --help` or see [CLI Reference](../reference/cli_reference.md).

## What's Next?

Now that you have run your first analysis, here is a suggested progression through the documentation:

- **[Running the Workflow](../usage/running_workflow.md)** -- Advanced execution options, resuming failed runs, and cluster configuration
- **[Input Preparation](../usage/input_preparation.md)** -- Detailed file format specifications and tips for preparing your data
- **[Interpreting Results](../usage/interpreting_results.md)** -- In-depth guide to understanding scores, p-values, and outlier detection
- **[Visualization Guide](../usage/visualization_guide.md)** -- Plot allele frequency trajectories and generate publication-ready figures
- **[Configuration Reference](../reference/configuration.md)** -- Full documentation of every configuration option
- **[Tutorial](../examples/tutorial.md)** -- End-to-end walkthrough with example data
