Quickstart
==========

This quickstart guide will help you run AlleleFlux on your data.

Basic Workflow
-------------

AlleleFlux is designed to be run using Snakemake, which manages the workflow and dependencies between tasks. The workflow consists of two main steps:

1. **Step 1**: Profile samples and generate an eligibility table
2. **Step 2**: Analyze allele frequencies and calculate scores

Running the Workflow
----------------------

First, copy the configuration file template and modify it for your data:

.. code-block:: bash

    cd /path/to/AlleleFlux
    cp smk_workflow/config.yml my_config.yml

Edit ``my_config.yml`` to specify your input files, output directories, and analysis parameters.

Run Step 1 (profiling samples):

.. code-block:: bash

    cd smk_workflow
    snakemake -s step1.smk --configfile ../my_config.yml --profile cornell_profile/

Run Step 2 (analyzing alleles):

.. code-block:: bash

    snakemake -s step2.smk --configfile ../my_config.yml --profile cornell_profile/

Command-line Tools
-----------------

AlleleFlux also provides command-line tools that you can use directly:

Profiling MAGs:

.. code-block:: bash

    alleleflux-profile --bam_path /path/to/bam --fasta_path /path/to/fasta \
    --prodigal_fasta /path/to/genes.fna --output_dir /path/to/output

Analyzing allele frequencies:

.. code-block:: bash

    alleleflux-allele-freq --magID MAG_ID --mag_metadata_file metadata.tsv \
    --fasta reference.fa --breath_threshold 0.1 --data_type single \
    --output_dir /path/to/output

For more detailed information on using AlleleFlux, please refer to the :doc:`../usage/running_workflow` guide.