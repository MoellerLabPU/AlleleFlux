# Installation

## System Requirements

Before installing AlleleFlux, ensure the following are available on your system:

| Requirement | Version | Notes |
|-------------|---------|-------|
| **Python** | >= 3.9 | Required by AlleleFlux core |
| **R** | >= 4.0 (recommended) | Required for CMH tests via rpy2; R packages: `stats` (built-in), `tidyr` |
| **samtools** | >= 1.10 (recommended) | For BAM file indexing and processing (provides `htslib` for pysam) |
| **Snakemake** | >= 8.0.0 | Workflow orchestration engine |

:::{note}
If you install via conda/mamba or `environment.yml`, all of these dependencies (including R and R packages) are handled automatically. Manual dependency management is only needed for pip-based installations.
:::

## Prerequisites

- **Python >= 3.9**
- **R** (required for CMH tests via rpy2; `r-base` and `r-tidyr` are included in the conda environment)
- A Unix-like operating system (Linux or macOS); Windows is supported for Python scripts but not for the Snakemake workflow

## Using Conda or Mamba (Recommended)

AlleleFlux is available on [Bioconda](https://anaconda.org/bioconda/alleleflux). Using conda (or mamba) ensures all dependencies -- including R and system libraries -- are installed correctly.

```bash
# Create a new environment with AlleleFlux
conda create -n alleleflux -c conda-forge -c bioconda alleleflux

# Activate the environment
conda activate alleleflux
```

Or install into an existing environment:

```bash
conda install -c conda-forge -c bioconda alleleflux
```

:::{tip}
We recommend [mamba](https://mamba.readthedocs.io/) as a faster drop-in replacement for conda:

```bash
mamba create -n alleleflux -c conda-forge -c bioconda alleleflux
```

:::

## Using environment.yml

If you want to install from the repository with all dependencies pre-configured (including R and R packages), use the provided `environment.yml`:

```bash
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux
conda env create -f environment.yml
conda activate alleleflux
```

The `environment.yml` file installs all dependencies from `conda-forge` and `bioconda`, including:

- All Python packages (pysam, rpy2, snakemake, etc.)
- R base and the `tidyr` R package (required for CMH tests)
- AlleleFlux itself in editable/development mode via `pip install -e .`

This is the recommended approach if you plan to modify the source code or need a reproducible environment that matches the development setup exactly.

## From Source

For development or to use the latest features:

```bash
# Clone the repository
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux

# Option 1: Use the provided environment file (recommended)
conda env create -f environment.yml
conda activate alleleflux

# Option 2: Install with pip (requires R and system dependencies separately)
pip install -e .
```

The `environment.yml` file includes all Python and R dependencies needed by AlleleFlux.

## From PyPI

```bash
pip install alleleflux
```

:::{warning}
Installing via pip does **not** install R or system-level dependencies (e.g., `htslib` for pysam). You must ensure R, r-tidyr, and samtools/htslib are available in your environment. For most users, the conda installation is strongly recommended.
:::

## Verify Installation

Confirm AlleleFlux is installed correctly:

```bash
# Check version
alleleflux --version

# Show installation info
alleleflux info

# Verify a CLI tool is available
alleleflux-profile --help
```

Expected output from `alleleflux info`:

```text
AlleleFlux v0.1.4
Package location: /path/to/alleleflux
Snakefile: /path/to/alleleflux/smk_workflow/alleleflux_pipeline/Snakefile
```

## Dependencies

AlleleFlux depends on the following major packages (installed automatically):

| Package | Purpose |
|---------|---------|
| `snakemake` >= 8.0 | Workflow orchestration |
| `pysam` | BAM file processing |
| `pandas`, `numpy` | Data manipulation |
| `scipy`, `statsmodels` | Statistical tests and LMM |
| `rpy2` | R integration for CMH tests |
| `biopython` | Sequence parsing (Prodigal genes) |
| `matplotlib`, `seaborn` | Visualization |
| `tqdm` | Progress bars |
| `click`, `questionary` | CLI interface |

For the complete list, see [`pyproject.toml`](https://github.com/MoellerLabPU/AlleleFlux/blob/main/pyproject.toml).

## Troubleshooting Installation

### rpy2 / R Issues

If you see errors related to `rpy2` or R not being found:

```text
rpy2.rinterface_lib.openrlib.ffi.error: R_HOME not found
```

**Solution**: Install R via conda rather than relying on a system R installation:

```bash
conda install -c conda-forge r-base r-tidyr
```

Ensure `R_HOME` is set correctly. Within a conda environment, this is handled automatically. If installing outside conda, set it manually:

```bash
export R_HOME=$(R RHOME)
```

### pysam Compilation Issues

If `pysam` fails to compile during `pip install` (missing `htslib` headers):

```text
fatal error: htslib/sam.h: No such file or directory
```

**Solution**: Install `pysam` from conda-forge instead of building from source:

```bash
conda install -c conda-forge pysam
```

Alternatively, install `htslib` and `samtools` system-wide and ensure they are on your `PATH` before running `pip install`.

### Snakemake Version Conflicts

AlleleFlux requires Snakemake >= 8.0.0. If you have an older version installed:

```bash
# Check current version
snakemake --version

# Upgrade within conda
conda install -c conda-forge -c bioconda 'snakemake>=8.0.0'
```

Snakemake 8.x introduced a new plugin-based executor architecture. AlleleFlux requires the `snakemake-executor-plugin-cluster-generic` package for cluster execution, which is included in the conda and `environment.yml` installations.

## Available CLI Tools

AlleleFlux provides 30+ command-line tools organized by function:

| Category | Tools |
|----------|-------|
| **Main** | `alleleflux run`, `alleleflux init`, `alleleflux info`, `alleleflux tools` |
| **Analysis** | `alleleflux-profile`, `alleleflux-allele-freq`, `alleleflux-scores`, `alleleflux-outliers` |
| **Preprocessing** | `alleleflux-qc`, `alleleflux-eligibility`, `alleleflux-metadata` |
| **Statistics** | `alleleflux-lmm`, `alleleflux-cmh`, `alleleflux-two-sample-paired`, `alleleflux-single-sample` |
| **Evolution** | `alleleflux-dnds-from-timepoints` |
| **Visualization** | `alleleflux-plot-trajectories`, `alleleflux-track-alleles`, `alleleflux-terminal-nucleotide` |
| **Accessory** | `alleleflux-create-mag-mapping`, `alleleflux-add-bam-path`, `alleleflux-list-mags` |

To list all available tools:

```bash
alleleflux tools
```

For a complete reference, see [CLI Reference](../reference/cli_reference.md).

## Next Steps

- [Quickstart](quickstart.md) -- Run your first analysis in minutes
- [Overview](overview.md) -- Understand the AlleleFlux pipeline
- [Input Preparation](../usage/input_preparation.md) -- Prepare your data files
