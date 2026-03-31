# AlleleFlux

[![Bioconda](https://img.shields.io/conda/vn/bioconda/alleleflux?label=bioconda)](https://anaconda.org/bioconda/alleleflux)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Documentation Status](https://readthedocs.org/projects/alleleflux/badge/?version=latest)](https://alleleflux.readthedocs.io/en/latest/?badge=latest)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/MoellerLabPU/AlleleFlux)

**AlleleFlux** is a bioinformatics toolkit for analyzing allele frequency changes in metagenomic time-series data. While tools exist for tracking strain-level variation in microbial communities, there has been a gap in methods for systematically detecting natural selection from longitudinal metagenomics. AlleleFlux fills this gap by providing a reproducible, end-to-end pipeline that quantifies parallel and divergent evolutionary dynamics across microbial populations, enabling researchers to pinpoint genomic targets of selection directly from shotgun metagenomic data.

## Features

- **Single-timepoint and longitudinal study designs** -- flexible analysis for cross-sectional or time-series experiments
- **Parallelism and divergence scoring** -- quantify parallel allele frequency changes across replicates and divergence between groups, with outlier gene detection
- **Five statistical test types** -- two-sample paired/unpaired, single-sample, Linear Mixed Models (LMM), and Cochran-Mantel-Haenszel (CMH) tests
- **dN/dS ratio calculation** -- measure selection pressure on genes via the Nei-Gojobori method with path averaging
- **Automated quality control** -- coverage-based filtering and MAG eligibility determination
- **Taxonomic aggregation** -- aggregate results from phylum to species using GTDB taxonomy
- **Visualization workflow** -- generate allele frequency trajectory plots
- **30+ standalone CLI tools** -- run individual analysis steps outside the workflow
- **Snakemake-based workflow** -- reproducible, parallelized execution with SLURM/HPC support
- **Bioconda and PyPI distribution** -- install via conda, pip, or from source

These scores enable direct comparisons of evolutionary dynamics across taxa, genomes, and genes, helping identify loci under strong selection.

## Requirements

- **Python** >= 3.9
- **R** (required for CMH tests via `rpy2`)
- **samtools** (for BAM indexing)

All dependencies are bundled in the provided `environment.yml` for reproducible setup.

## Installation

### From Bioconda (Recommended)

```bash
# Install with conda (or mamba)
conda install -c conda-forge -c bioconda alleleflux

# Activate the environment
conda activate alleleflux
```

### Using environment.yml (Recommended for New Users)

```bash
# Create environment with all dependencies
conda env create -f environment.yml
conda activate alleleflux
```

### From PyPI

```bash
pip install alleleflux
```

### From Source

```bash
# Clone the repository
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux

# Create environment with dependencies
conda env create -f environment.yml
conda activate alleleflux

# Or install directly with pip
pip install -e .
```

## Input Files

| File | Description |
|------|-------------|
| Reference FASTA | Combined MAG contigs (header format: `<MAG_ID>.fa_<contig_ID>`) |
| Prodigal genes | Nucleotide ORF predictions (`.fna`) matching reference contig IDs |
| MAG mapping | TSV with columns: `mag_id`, `contig_id` |
| Metadata TSV | Sample info with columns: `sample_id`, `bam_path`, `group`, `time` |
| GTDB taxonomy | `gtdbtk.bac120.summary.tsv` from GTDB-Tk |

See [Input Preparation Guide](https://alleleflux.readthedocs.io/en/latest/usage/input_preparation.html) for detailed format specifications.

## Quick Start

### 1. Initialize Configuration

```bash
# Interactive configuration wizard
alleleflux init

# Or copy the template manually
cp $(python -c "import alleleflux; print(alleleflux.__path__[0])")/smk_workflow/config.template.yml config.yml
```

### 2. Edit Configuration

Edit `config.yml` with your paths and parameters. Here is the complete configuration with all options:

```yaml
run_name: "alleleflux_analysis"

# Input Files
input:
  fasta_path: ""                  # Reference FASTA file (required)
  prodigal_path: ""               # Prodigal nucleic acid output (.fna)
  metadata_path: ""               # Sample metadata file
  gtdb_path: ""                   # GTDB taxonomy file
  mag_mapping_path: ""            # MAG-to-contig mapping file

# Output Directory
output:
  root_dir: "./alleleflux_output"

log_level: "INFO"                 # DEBUG, INFO, WARNING, ERROR

# Analysis Configuration
analysis:
  data_type: "longitudinal"       # "single" or "longitudinal"
  allele_analysis_only: False     # Skip scoring/outlier detection
  use_lmm: True                   # Linear Mixed Models
  use_significance_tests: True    # Two-sample and single-sample tests
  use_cmh: True                   # CMH test

  timepoints_combinations:
    - timepoint: ["pre", "post"]
      focus: "post"

  groups_combinations:
    - ["treatment", "control"]

# Quality Control
quality_control:
  min_sample_num: 4
  breadth_threshold: 0.1
  coverage_threshold: 1
  disable_zero_diff_filtering: False

# Profiling
profiling:
  ignore_orphans: False
  min_base_quality: 30
  min_mapping_quality: 2
  ignore_overlaps: True

# Statistics
statistics:
  filter_type: "t-test"           # "t-test", "wilcoxon", or "both"
  preprocess_between_groups: True
  preprocess_within_groups: True
  max_zero_count: 4
  p_value_threshold: 0.05
  fdr_group_by_mag_id: False
  min_positions_after_preprocess: 1

# dN/dS Analysis
dnds:
  p_value_column: "q_value"
  dn_ds_test_type: "two_sample_paired_tTest"

# Compute Resources
resources:
  threads_per_job: 16
  mem_per_job: "8G"
  time: "24:00:00"
```

#### Configuration Parameters

| Section | Parameter | Description |
|---------|-----------|-------------|
| **input** | `fasta_path` | Reference FASTA with combined MAG contigs |
| | `prodigal_path` | Prodigal nucleotide predictions (`.fna`) |
| | `metadata_path` | Sample metadata TSV |
| | `gtdb_path` | GTDB-Tk taxonomy file |
| | `mag_mapping_path` | MAG-to-contig mapping TSV |
| **analysis** | `data_type` | `"longitudinal"` (multiple timepoints) or `"single"` |
| | `allele_analysis_only` | Skip significance tests, scoring, and outlier detection if `True` |
| | `use_lmm` | Enable Linear Mixed Models |
| | `use_significance_tests` | Enable two-sample/single-sample tests |
| | `use_cmh` | Enable Cochran-Mantel-Haenszel tests |
| | `timepoints_combinations` | Timepoint pairs with focus timepoint |
| | `groups_combinations` | Groups to compare |
| **quality_control** | `min_sample_num` | Minimum samples required per MAG |
| | `breadth_threshold` | Minimum coverage breadth (0-1) |
| | `coverage_threshold` | Minimum average coverage depth |
| **profiling** | `min_base_quality` | Minimum base quality score |
| | `min_mapping_quality` | Minimum mapping quality score |
| **statistics** | `filter_type` | Preprocessing filter type |
| | `p_value_threshold` | Significance threshold |
| | `fdr_group_by_mag_id` | Apply FDR correction per MAG |
| **dnds** | `p_value_column` | `"min_p_value"` or `"q_value"` |
| | `dn_ds_test_type` | Test type for filtering dN/dS results |
| **resources** | `threads_per_job` | Threads allocated per job |
| | `mem_per_job` | Memory per job (e.g., `"8G"`, `"100G"`) |
| | `time` | Wall time limit (HH:MM:SS) |

See [Configuration Reference](https://alleleflux.readthedocs.io/en/latest/reference/configuration.html) for complete documentation.

### 3. Run the Pipeline

```bash
# Run locally
alleleflux run --config config.yml --threads 16

# Dry run to preview jobs
alleleflux run --config config.yml --dry-run
```

#### Running on SLURM

For HPC clusters, copy the SLURM profile from the source repository:

```bash
# Copy SLURM profile (if installed from source)
cp -r $(python -c "import alleleflux; print(alleleflux.__path__[0])")/smk_workflow/slurm_profile ./

# Run with SLURM
alleleflux run --config config.yml --profile ./slurm_profile
```

The SLURM profile automatically submits jobs via `sbatch` with resources from your config.

## How It Works

AlleleFlux is powered by a **Snakemake workflow** that orchestrates the complete analysis:

```text
Input Files              Profile & QC           Statistical Analysis
━━━━━━━━━━━              ━━━━━━━━━━━━           ━━━━━━━━━━━━━━━━━━━━
• Reference FASTA        • Extract alleles      • Two-sample tests
• Prodigal genes         • Quality control      • LMM / CMH tests
• Metadata TSV           • Eligibility checks   • dN/dS calculation
• MAG mapping                    ↓
                         ┌─────────────────┐
                         │ Scoring & Viz   │
                         ├─────────────────┤
                         │ • Parallelism   │
                         │ • Divergence    │
                         │ • Outliers      │
                         │ • Trajectories  │
                         └─────────────────┘
```

**Pipeline Steps:**

1. **Profiling** -- Extract allele frequencies from BAM files for each MAG
2. **Quality Control** -- Filter samples by coverage breadth; determine MAG eligibility
3. **Statistical Testing** -- Apply appropriate tests based on experimental design
4. **Scoring** -- Calculate parallelism/divergence scores and identify outlier genes
5. **dN/dS Analysis** -- Calculate evolutionary rates for genes under selection

The workflow:

- Automatically parallelizes across samples and MAGs
- Handles checkpointing and restarts gracefully
- Supports local execution and HPC clusters (SLURM)
- Tracks provenance and ensures reproducibility

## Output

Results are organized in the output directory:

```text
alleleflux_output/
├── profiles/           # Per-sample allele frequency profiles
├── metadata/           # Per-MAG metadata tables
├── eligibility/        # MAG eligibility tables
├── allele_analysis/    # Allele frequency analysis results
├── significance_tests/ # Statistical test results (LMM, CMH, t-tests)
├── scores/             # Parallelism and divergence scores
├── outliers/           # Genes with high scores (selection targets)
└── dnds/               # dN/dS analysis results
```

See [Output Reference](https://alleleflux.readthedocs.io/en/latest/reference/outputs.html) for file format details.

## CLI Tools

AlleleFlux provides 30+ standalone command-line tools:

```bash
# List all available tools
alleleflux tools

# Main commands
alleleflux run --help          # Run the full pipeline
alleleflux init --help         # Interactive configuration
alleleflux info                # Show installation info

# Individual analysis tools
alleleflux-profile --help      # Profile MAGs from BAM files
alleleflux-qc --help           # Quality control
alleleflux-scores --help       # Calculate parallelism/divergence scores
alleleflux-dnds-from-timepoints --help  # Calculate dN/dS ratios
```

See [CLI Reference](https://alleleflux.readthedocs.io/en/latest/reference/cli_reference.html) for the complete list.

## Citing AlleleFlux

If you use AlleleFlux in your research, please cite:

> Uppal, S. & Moeller, A.H. AlleleFlux: A bioinformatics toolkit for analyzing allele frequency dynamics in metagenomic data. GitHub: https://github.com/MoellerLabPU/AlleleFlux

## Documentation

Full documentation: **[alleleflux.readthedocs.io](https://alleleflux.readthedocs.io/)**

- [Installation Guide](https://alleleflux.readthedocs.io/en/latest/getting_started/installation.html)
- [Quickstart Tutorial](https://alleleflux.readthedocs.io/en/latest/getting_started/quickstart.html)
- [Configuration Reference](https://alleleflux.readthedocs.io/en/latest/reference/configuration.html)
- [Interpreting Results](https://alleleflux.readthedocs.io/en/latest/usage/interpreting_results.html)

## Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

AlleleFlux is licensed under the [GNU General Public License v3.0](LICENSE).

## Acknowledgments

AlleleFlux was developed at [Princeton University](https://www.princeton.edu/) by Siddhartha Uppal and Andrew Moeller in the [Moeller Lab](https://github.com/MoellerLabPU).
