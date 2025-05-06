Input Preparation
==================

AlleleFlux requires several input files to run. This guide explains how to prepare each required input.

Required Inputs
----------------

1. **BAM Files**
   
   Sorted and indexed BAM files from metagenomic samples aligned to reference MAGs (Metagenome-Assembled Genomes).
   
   * File format: ``.bam`` (with accompanying ``.bam.bai`` index file)
   * Each BAM file should contain alignments of reads to reference MAGs
   * BAM files can have any naming convention, as they are referenced through the metadata file

2. **Reference FASTA**
   
   A FASTA file containing the reference MAG sequences.
   
   * File format: ``.fa`` or ``.fasta``
   * Should contain all contigs referenced in the BAM files

3. **Prodigal Genes FASTA**
   
   A FASTA file containing predicted genes from Prodigal.
   
   * File format: ``.fna``
   * Should be generated using Prodigal on the reference FASTA

4. **Metadata File**
   
   A tab-separated file containing metadata for each sample, including BAM file paths.
   
   * File format: ``.txt`` or ``.tsv``
   * Required columns:
     - ``sample_id``: Unique identifier for each sample
     - ``bam_path``: Full path to the BAM file for the sample
     - ``subjectID``: Identifier for the subject/replicate
     - ``group``: Group identifier (e.g., "control", "treatment")
     - ``replicate``: Replicate identifier
   * For longitudinal data, also include:
     - ``time``: Timepoint identifier (e.g., "pre", "post", "day1")

5. **GTDB Taxonomy File** (Optional)
   
   A taxonomy file from GTDB (Genome Taxonomy Database) for taxonomic classification.
   
   * File format: ``.tsv``
   * Used for taxonomic analysis

Example Metadata File
---------------------

.. code-block:: text

    sample_id    bam_path                               subjectID    group      replicate    time
    S1           /path/to/bamfiles/S1.bam               mouse1       control    A           pre
    S2           /path/to/bamfiles/sample2.bam          mouse2       control    B           pre
    S3           /path/to/bamfiles/patient3_t1.bam      mouse3       treatment  A           pre
    S4           /path/to/bamfiles/S4_rep2.bam          mouse4       treatment  B           pre
    S5           /path/to/bamfiles/S5.sorted.bam        mouse1       control    A           post
    S6           /path/to/bamfiles/S6_aligned.bam       mouse2       control    B           post
    S7           /path/to/bamfiles/S7_final.bam         mouse3       treatment  A           post
    S8           /path/to/bamfiles/S8_processed.bam     mouse4       treatment  B           post

Configuration File
-------------------

Update the ``config.yml`` file with the paths to your input files:

.. code-block:: yaml

    # Inputs
    input:
      bam_dir: "/path/to/bamfiles"  # Used for backward compatibility
      fasta_path: "/path/to/reference.fa"
      prodigal_path: "/path/to/prodigal_genes.fna"
      metadata_path: "/path/to/metadata.txt"  # Must include bam_path column
      gtdb_path: "/path/to/gtdb_taxonomy.tsv"
    
    # Outputs
    output:
      root_dir: "/path/to/output_directory"
    
    # Parameters
    timepoints_combinations:
      - ["pre", "post"]
    
    groups_combinations:
      - ["control", "treatment"]
    
    # Analysis options
    analysis_options:
      use_lmm: True
      use_significance_tests: True
      use_cmh: True                # Enable Cochran-Mantel-Haenszel tests
      data_type: "longitudinal"  # or "single"
    
    # Focus timepoints for CMH test
    focus_timepoints:
      pre_post: "post"  # Specify which timepoint to focus on for CMH analysis
    
    min_sample_num: 4
    breadth_threshold: 0.1
    disable_zero_diff_filtering: False
    alpha: 0.05
    test_type: "both"
    preprocess_two_sample: True
    alpha: 0.05
    test_type: "both"
    preprocess_two_sample: True

Adding BAM Paths to Metadata
---------------------------

AlleleFlux includes a utility script to help you add BAM file paths to your existing metadata file:

.. code-block:: bash

    python alleleflux/accessory/add_bam_path_to_metadata.py \
        --metadata /path/to/metadata.tsv \
        --output /path/to/updated_metadata.tsv \
        --bam-dir /path/to/bamfiles \
        --bam-extension .bam

This script will:

1. Read your existing metadata file
2. Search for BAM files in the specified directory that match each sample ID
3. Add a ``bam_path`` column to your metadata file
4. Save the updated metadata to a new file

You can then use this updated metadata file with AlleleFlux.

Options:
  * ``--metadata``: Path to your existing metadata file (required)
  * ``--output``: Path to save the updated metadata file (required)
  * ``--bam-dir``: Directory containing BAM files (default: current directory)
  * ``--bam-extension``: Extension of BAM files (default: .bam)
  * ``--drop-missing``: Drop samples without matching BAM files (optional)

Focus Timepoints for CMH Test
----------------------------

The Cochran-Mantel-Haenszel (CMH) test requires a focus timepoint to be specified for each timepoint combination. The focus timepoint is used to identify significant allele frequency changes between groups at a specific timepoint.

In the configuration file, specify the focus timepoint for each timepoint combination:
The CMH test will calculate significance scores based on the specified focus timepoint. This is particularly useful for identifying evolutionary changes at specific points in time.