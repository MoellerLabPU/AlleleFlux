# AlleleFlux Example Data

This directory contains minimal synthetic data for testing and learning AlleleFlux.

## Contents

```
example_data/
├── README.md                    # This file
├── reference/
│   ├── combined_mags.fasta      # Reference FASTA with 2 MAGs
│   ├── prodigal_genes.fna       # Gene predictions
│   ├── mag_mapping.tsv          # Contig-to-MAG mapping
│   └── gtdbtk_taxonomy.tsv      # Mock GTDB taxonomy
├── metadata/
│   └── sample_metadata.tsv      # Sample metadata (8 samples)
├── profiles/                    # Pre-generated profile files
│   ├── control_subj1_pre/
│   │   ├── control_subj1_pre_TEST_MAG_001_profiled.tsv.gz
│   │   └── control_subj1_pre_TEST_MAG_002_profiled.tsv.gz
│   └── ... (8 sample directories)
├── significant_sites/
│   └── significant_sites.tsv    # Example significant sites for visualization
└── config_example.yml           # Working configuration file
```

## Dataset Overview

- **2 MAGs**: `TEST_MAG_001` and `TEST_MAG_002`
- **8 samples**: 4 control, 4 treatment (2 timepoints each: pre/post)
- **4 subjects**: Each with samples at both timepoints (longitudinal design)
- **~2000 positions per MAG**: Sufficient to demonstrate allele frequency analysis

## Usage

### Running the Full Pipeline

```bash
cd /path/to/AlleleFlux
alleleflux run --config docs/source/examples/example_data/config_example.yml
```

### Running Individual Steps

```bash
# Profile a sample (already done - profiles provided)
alleleflux-profile \
    --bam_path /path/to/your.bam \
    --fasta_path docs/source/examples/example_data/reference/combined_mags.fasta \
    --prodigal_fasta docs/source/examples/example_data/reference/prodigal_genes.fna \
    --mag_mapping_file docs/source/examples/example_data/reference/mag_mapping.tsv \
    --output_dir output/profiles

# Run visualization workflow
alleleflux-terminal-nuc-analysis \
    --significant_sites docs/source/examples/example_data/significant_sites/significant_sites.tsv \
    --profile_dir docs/source/examples/example_data/profiles \
    --metadata docs/source/examples/example_data/metadata/sample_metadata.tsv \
    --group treatment \
    --timepoint post \
    --output results/terminal
```

## Generating Larger Datasets

For testing with larger datasets, use the provided generation script:

```bash
python docs/source/examples/generate_synthetic_data.py \
    --num_mags 10 \
    --num_samples 20 \
    --num_positions 5000 \
    --output_dir my_test_data
```

See `generate_synthetic_data.py --help` for all options.

## File Format Specifications

### Profile Files (`*_profiled.tsv.gz`)

| Column | Type | Description |
|--------|------|-------------|
| `contig` | string | Contig identifier |
| `position` | int | 0-based genomic position |
| `gene_id` | string | Overlapping gene ID |
| `ref_base` | string | Reference base (A/C/G/T) |
| `A` | int | Count of adenine bases |
| `T` | int | Count of thymine bases |
| `G` | int | Count of guanine bases |
| `C` | int | Count of cytosine bases |

### Metadata File (`sample_metadata.tsv`)

| Column | Description |
|--------|-------------|
| `sample_id` | Unique sample identifier |
| `bam_path` | Path to BAM file (or profile directory) |
| `subjectID` | Biological replicate/subject ID |
| `group` | Experimental group (control/treatment) |
| `replicate` | Replicate letter within group |
| `time` | Timepoint (pre/post) |

### Significant Sites (`significant_sites.tsv`)

| Column | Description |
|--------|-------------|
| `mag_id` | MAG identifier |
| `contig` | Contig identifier |
| `position` | 0-based position |
| `gene_id` | Gene identifier |
| `test_type` | Statistical test used |
| `min_p_value` | Minimum p-value across tests |
| `q_value` | FDR-adjusted p-value |
