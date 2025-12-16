# Installation

## Prerequisites

AlleleFlux requires Python 3.7+ and depends on standard bioinformatics packages (pandas, numpy, scipy, biopython, pysam, rpy2, etc.). All dependencies are managed automatically during installation.

## Using Conda or Mamba (Recommended)

You can use either `conda` or `mamba` (both are fine). AlleleFlux is available on conda-forge and Bioconda:

```bash
# Install with conda (or mamba)
conda install -c conda-forge -c bioconda alleleflux

# Activate the environment
conda activate alleleflux
```

## From Source

For development or if you prefer installing from source:

```bash
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux
pip install -e .
```

This installs AlleleFlux in editable mode with all command-line tools available.

## Verify Installation

Confirm AlleleFlux is installed correctly:

```bash
alleleflux --help
alleleflux-profile --help
```

You should see help text describing available commands and options.

## Available CLI Tools

AlleleFlux provides 30+ command-line tools organized by function:

- **Analysis**: `alleleflux-profile`, `alleleflux-allele-freq`, `alleleflux-scores`, `alleleflux-outliers`, `alleleflux-dnds`
- **Preprocessing**: `alleleflux-qc`, `alleleflux-eligibility`, `alleleflux-metadata`
- **Statistics**: `alleleflux-lmm`, `alleleflux-cmh`, `alleleflux-two-sample-paired`, `alleleflux-single-sample`
- **Visualization**: `alleleflux-prepare-metadata`, `alleleflux-terminal-nucleotide`, `alleleflux-track-alleles`, `alleleflux-plot-trajectories`
- **Utilities**: `alleleflux-create-mag-mapping`, `alleleflux-add-bam-path`

For a complete reference, see [CLI Reference](../reference/cli_reference.md).

## Next Steps

Proceed to [Quickstart](quickstart.md) to run your first analysis, or read the [Overview](overview.md) for a conceptual introduction.
