
[![Documentation Status](https://readthedocs.org/projects/alleleflux/badge/?version=latest)](https://alleleflux.readthedocs.io/en/latest/?badge=latest) [![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/MoellerLabPU/AlleleFlux)


# AlleleFlux

A toolkit for analyzing allele frequency changes and dN/dS ratios in metagenomic data.

## Installation

You can use either `conda` or `mamba` (both are fine).

```bash
# Install with conda (or mamba)
conda install -c conda-forge -c bioconda alleleflux

# Activate the environment (if a new env was created by your conda setup)
conda activate alleleflux
```

### From source

If you prefer to install from source:

```bash
git clone https://github.com/MoellerLabPU/AlleleFlux.git
cd AlleleFlux
pip install -e .
```

## Usage

### Run the Pipeline

```bash
# Configure your analysis
cp alleleflux/smk_workflow/config.template.yml config.yml
# Edit config.yml with your input paths and parameters

# Run
alleleflux run --config config.yml --threads 16
```

### Input Files

| File | Description |
|------|-------------|
| BAM files | Sorted and indexed (`*.sorted.bam` + `*.sorted.bam.bai`) |
| Reference FASTA | Combined MAG contigs (header: `<MAG_ID>.fa_<contig_ID>`) |
| Prodigal genes | Nucleotide ORF predictions matching reference IDs |
| Metadata TSV | Columns: `sample_id`, `bam_path`, `replicate`, `subjectID`, `time`, `group` |
| GTDB taxonomy | Optional: `gtdbtk.bac120.summary.tsv` from GTDB-Tk |

### CLI Tools

```bash
# See all available commands
alleleflux --help

# Examples
alleleflux-profile --help      # Profile MAGs from BAM files
alleleflux-qc --help           # Quality control
alleleflux-scores --help       # Calculate parallelism/divergence scores
alleleflux-dnds --help         # Calculate dN/dS ratios from timepoint comparisons
```

## Documentation

Full documentation: [alleleflux.readthedocs.io](https://alleleflux.readthedocs.io/)

## License

See [LICENSE](LICENSE) for details.
