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

Edit ``my_config.yml`` to specify your input files, output directories, and analysis parameters. Key configuration options include:

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
    
    # Analysis settings
    analysis:
      use_significance_tests: true  # Enable two-sample and single-sample tests
      use_lmm: true                # Enable Linear Mixed Models
      use_cmh: true                # Enable Cochran-Mantel-Haenszel tests

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

**Preprocessing utilities:**

Create MAG mapping files (run before main workflow):

.. code-block:: bash

    alleleflux-create-mag-mapping --dir /path/to/mag_fastas \
    --extension fa --output-fasta combined.fasta \
    --output-mapping mapping.tsv

Add BAM paths to metadata:

.. code-block:: bash

    alleleflux-add-bam-path --metadata metadata.tsv \
    --bam-dir /path/to/bams --output updated_metadata.tsv

**Profiling MAGs:**

.. code-block:: bash

    alleleflux-profile --bam_path /path/to/bam --fasta_path /path/to/fasta \
    --prodigal_fasta /path/to/genes.fna --output_dir /path/to/output

**Running statistical tests:**

CMH Test:

.. code-block:: bash

    alleleflux-cmh --input_df /path/to/longitudinal.tsv.gz \
    --preprocessed_df /path/to/preprocessed.tsv.gz \
    --min_sample_num 4 --mag_id MAG_ID --data_type longitudinal \
    --cpus 16 --output_dir /path/to/output

Linear Mixed Models:

.. code-block:: bash

    alleleflux-lmm --input_df /path/to/data.tsv.gz \
    --preprocessed_df /path/to/preprocessed.tsv.gz \
    --group group_name --mag_id MAG_ID --cpus 16 \
    --output_dir /path/to/output

**Analyzing allele frequencies:**

.. code-block:: bash

    alleleflux-allele-freq --magID MAG_ID --mag_metadata_file metadata.tsv \
    --fasta reference.fa --breadth_threshold 0.1 --data_type single \
    --output_dir /path/to/output

**Quality control:**

.. code-block:: bash

    alleleflux-qc --rootDir /path/to/profiles --fasta reference.fa \
    --breadth_threshold 0.1 --output_dir /path/to/qc

**Eligibility table generation:**

.. code-block:: bash

    alleleflux-eligibility --qc_dir /path/to/qc --min_sample_num 4 \
    --output_file eligibility.tsv --data_type longitudinal

For help with any command, use the ``-h`` or ``--help`` flag:

.. code-block:: bash

    alleleflux-profile --help
    alleleflux-cmh --help

For more detailed information on using AlleleFlux, please refer to the :doc:`../usage/running_workflow` guide.