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
2. Perform statistical tests based on your configuration
3. Calculate parallelism and divergence scores
4. Identify genes with exceptionally high scores (outliers)

Customizing the Workflow
-------------------------

You can customize the workflow by editing the ``config.yml`` file (see :doc:`input_preparation` for details).

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

Checkpoint Files
~~~~~~~~~~~~~~~~

The workflow creates checkpoint files at various stages. If you need to restart a failed run, Snakemake will automatically pick up from the last successful checkpoint.

Troubleshooting
---------------

If you encounter issues when running the workflow:

1. Check the Snakemake log files in the ``smk_workflow/logs/`` directory
2. Ensure that all input files are in the correct format
3. Verify that you have sufficient resources (memory, CPU, disk space)
4. Check that all dependencies are installed correctly