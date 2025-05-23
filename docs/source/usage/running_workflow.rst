Running the Workflow
=====================

AlleleFlux uses Snakemake to manage the workflow. This guide explains how to run the workflow and understand the different steps.

Workflow Overview
------------------

The AlleleFlux workflow consists of two main steps:

1. **Step 1**: Profile samples and generate an eligibility table
   - Profile BAM files to extract allele frequencies
   - Generate MAG metadata
   - Perform quality control
   - Create eligibility tables

2. **Step 2**: Analyze alleles and calculate scores
   - Analyze allele frequencies
   - Perform statistical tests
   - Calculate scores
   - Identify outlier genes

Running with Snakemake
------------------------

Prerequisites
~~~~~~~~~~~~~~

- AlleleFlux installed (see :doc:`../getting_started/installation`)
- Snakemake installed
- Input files prepared (see :doc:`input_preparation`)

Step 1: Profiling Samples
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the first step of the workflow:

.. code-block:: bash

    cd /path/to/AlleleFlux/smk_workflow
    snakemake -s step1.smk --profile cornell_profile/

This step will:

1. Process each BAM file to extract allele frequencies for each MAG
2. Generate metadata for each MAG based on your input metadata file
3. Perform quality control to filter out low-quality samples
4. Create eligibility tables to determine which MAGs can be used for statistical tests

Step 2: Analyzing Alleles
~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the second step of the workflow:

.. code-block:: bash

    cd /path/to/AlleleFlux/smk_workflow
    snakemake -s step2.smk --profile cornell_profile/

This step will:

1. Analyze allele frequencies for each eligible MAG
2. Perform statistical tests based on your configuration:
   - Two-sample tests (paired and unpaired) for comparing groups
   - Single-sample tests for analyzing within-group changes
   - Linear Mixed Models (LMM) for longitudinal data analysis
   - Cochran-Mantel-Haenszel (CMH) tests for stratified analysis of count data
3. Calculate parallelism and divergence scores
4. Identify genes with exceptionally high scores (outliers)

Available Statistical Tests
---------------------------

AlleleFlux supports several statistical testing approaches:

Two-Sample Tests
~~~~~~~~~~~~~~~~

These tests compare allele frequencies between groups:

- **Unpaired**: For comparing independent samples from different groups
- **Paired**: For comparing matched samples (e.g., before and after treatment)

Single-Sample Tests
~~~~~~~~~~~~~~~~~~

These tests analyze changes within a single group over time.

Linear Mixed Models (LMM)
~~~~~~~~~~~~~~~~~~~~~~~~

LMM tests account for complex experimental designs with fixed and random effects, particularly useful for longitudinal data.

Cochran-Mantel-Haenszel (CMH) Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CMH test is a stratified analysis of count data that:

- Tests for association between allele changes and conditions while controlling for confounding factors
- Provides position-by-position assessment of allele frequency changes
- Is especially powerful for detecting parallel evolutionary changes

To enable/disable specific tests, modify these settings in the ``config.yml`` file:

.. code-block:: yaml

    analysis:
      use_significance_tests: true  # Enable/disable two-sample and single-sample tests
      use_lmm: true                # Enable/disable Linear Mixed Models
      use_cmh: true                # Enable/disable Cochran-Mantel-Haenszel tests

Customizing the Workflow
-------------------------

You can customize the workflow by editing the ``config.yml`` file (see :doc:`input_preparation` for details). Key configuration options include:

.. code-block:: yaml

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

Advanced Usage
---------------

Running on a Compute Cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AlleleFlux supports running on a compute cluster through Snakemake's cluster support. To run on a cluster:

1. Create a cluster profile for your system (see Snakemake documentation)
2. Run the workflow using your cluster profile:

.. code-block:: bash

    snakemake -s step1.smk --profile your_cluster_profile/
    snakemake -s step2.smk --profile your_cluster_profile/

Example cluster profiles are provided for:

- Cornell BioHPC (``cornell_profile/``)
- Princeton Della (``della_profile/``)

You can adapt these profiles for your specific computing environment.

Output Files and Directories
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The workflow generates several output directories:

.. code-block:: text

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

Checkpoint Files
~~~~~~~~~~~~~~~~

The workflow creates checkpoint files at various stages. If you need to restart a failed run, Snakemake will automatically pick up from the last successful checkpoint.

Troubleshooting
---------------

If you encounter issues when running the workflow:

1. Check the Snakemake log files in the ``logs/`` directory
2. Ensure that all input files are in the correct format
3. Verify that you have sufficient resources (memory, CPU, disk space)
4. Check that all dependencies are installed correctly

Common Issues
~~~~~~~~~~~~~

- **Error in rpy2 or R dependencies**: Ensure you have R installed and R packages required for CMH tests (e.g., stats)
- **Memory errors**: Increase the memory allocation in your Snakemake profile
- **Missing files**: Check paths in your config.yml file