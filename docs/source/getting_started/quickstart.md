# Quickstart

Get started with AlleleFlux in minutes.

## Prerequisites

Ensure you have:

- AlleleFlux installed (see [Installation](installation.md))

- Input files prepared:

  - BAM files (sorted and indexed)
  - Reference FASTA (MAG contigs)
  - Prodigal gene predictions
  - Sample metadata TSV (with `bam_path` column)
  - MAG mapping file (contig → MAG assignments)

## Configuration

Copy the configuration template and edit for your data:

```bash
cp alleleflux/smk_workflow/config.template.yml my_config.yml
```

Minimal configuration:

```yaml
data_type: "longitudinal"  # or "single"

input:
  fasta_path: "reference.fa"
  prodigal_path: "genes.fna"
  metadata_path: "metadata.tsv"
  mag_mapping_path: "mag_mapping.tsv"

output:
  root_dir: "output/"

analysis:
  use_significance_tests: true
  use_lmm: true
  use_cmh: true
```

See [Configuration Reference](../reference/configuration.md) for all options.

## Run the Pipeline

```bash
alleleflux run --config my_config.yml --threads 16
```

This command profiles samples, runs quality control, performs statistical tests, calculates scores, and identifies outliers.

## Output

Results are organized in the output directory:

```text
output/
├── profiles/               # Allele frequency profiles
├── metadata/               # Per-MAG metadata
├── eligibility/            # MAG eligibility tables
├── allele_analysis/        # Frequency analysis results
├── significance_tests/     # Statistical test results
├── scores/                 # Parallelism/divergence scores
└── outliers/               # High-scoring genes
```

Key output files:

- `scores/processed/combined/scores_*.tsv` - MAG-level scores
- `outliers/*/outlier_genes_*.tsv` - Genes under selection
- `significance_tests/cmh/*.tsv.gz` - Position-level p-values

See [Output Files Reference](../reference/outputs.md) for complete specifications.

## Individual CLI Tools

Run individual analysis steps:

```bash
# Prepare MAG mapping
alleleflux-create-mag-mapping --dir mag_fastas/ --extension fa \
    --output-fasta combined.fasta --output-mapping mapping.tsv

# Profile a single sample
alleleflux-profile --bam_path sample.bam --fasta_path reference.fa \
    --prodigal_fasta genes.fna --output_dir profiles/

# Quality control
alleleflux-qc --rootDir profiles/ --fasta reference.fa \
    --breadth_threshold 0.1 --output_dir qc/
```

For complete CLI reference, see [CLI Reference](../reference/cli_reference.md).

## Next Steps

- [Running the Workflow](../usage/running_workflow.md) - Detailed workflow guide
- [Visualization Guide](../usage/visualization_guide.md) - Plot allele trajectories
- [Tutorial](../examples/tutorial.md) - Full tutorial with example data
- [Interpreting Results](../usage/interpreting_results.md) - Understanding your results
