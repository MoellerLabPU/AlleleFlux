# AlleleFlux Scripts

This directory contains scripts for the AlleleFlux pipeline, organized by functionality.

## Directory Structure

- **utilities/**: Core utility functions and shared modules
  - `utilities.py`: General utility functions used across multiple scripts
  - `supress_warning.py`: Warning suppression utilities

- **preprocessing/**: Scripts for data preparation and initial processing
  - `mag_metadata.py`: Handling MAG metadata
  - `quality_control.py`: Quality control procedures
  - `eligibility_table.py`: Generate eligibility tables
  - `preprocess_two_sample.py`: Preprocessing for two-sample analysis

- **analysis/**: Scripts for analyzing processed data
  - `profile_mags.py`: Profiling metagenome-assembled genomes
  - `allele_freq.py`: Analyzing allele frequencies
  - `gene_scores.py`: Computing gene scores
  - `taxa_scores.py`: Computing scores for taxonomic groups
  - `scores.py`: General scoring utilities
  - `outliers_genes.py`: Detection of outlier genes

- **statistics/**: Statistical testing and analysis scripts
  - `LMM.py`: Linear mixed model implementation
  - `single_sample.py`: Single sample statistical tests
  - `two_sample_paired.py`: Paired two-sample tests
  - `two_sample_unpaired.py`: Unpaired two-sample tests

## Usage

All scripts are referenced in the Snakemake workflow using relative paths in the config.yml file. 