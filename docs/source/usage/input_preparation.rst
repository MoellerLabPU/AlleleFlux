Input Preparation
==================

AlleleFlux requires several input files to run. This guide explains how to prepare each required input.

Required Inputs
----------------

1. **BAM Files**
   
   Sorted and indexed BAM files from metagenomic samples aligned to reference MAGs (Metagenome-Assembled Genomes).
   
   * File format: ``.sorted.bam`` (with accompanying ``.sorted.bam.bai`` index file)
   * Naming convention: ``{sample_id}.sorted.bam``
   * Each BAM file should contain alignments of reads to reference MAGs

2. **Reference FASTA**
   
   A FASTA file containing the reference MAG sequences.
   
   * File format: ``.fa`` or ``.fasta``
   * Should contain all contigs referenced in the BAM files

3. **Prodigal Genes FASTA**
   
   A FASTA file containing predicted genes from Prodigal.
   
   * File format: ``.fna``
   * Should be generated using Prodigal on the reference FASTA

4. **Metadata File**
   
   A tab-separated file containing metadata for each sample.
   
   * File format: ``.txt`` or ``.tsv``
   * Required columns:
     - ``sample_id``: Unique identifier for each sample
     - ``file_path``: Path to the BAM file for the sample
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

    sample_id    file_path                               subjectID    group      replicate    time
    S1           /path/to/bamfiles/S1.sorted.bam         mouse1       control    A           pre
    S2           /path/to/bamfiles/S2.sorted.bam         mouse2       control    B           pre
    S3           /path/to/bamfiles/S3.sorted.bam         mouse3       treatment  A           pre
    S4           /path/to/bamfiles/S4.sorted.bam         mouse4       treatment  B           pre
    S5           /path/to/bamfiles/S5.sorted.bam         mouse1       control    A           post
    S6           /path/to/bamfiles/S6.sorted.bam         mouse2       control    B           post
    S7           /path/to/bamfiles/S7.sorted.bam         mouse3       treatment  A           post
    S8           /path/to/bamfiles/S8.sorted.bam         mouse4       treatment  B           post

Configuration File
-------------------

Update the ``config.yml`` file with the paths to your input files:

.. code-block:: yaml

    # Inputs
    bamDir: "/path/to/bamfiles"
    fasta: "/path/to/reference.fa"
    prodigal: "/path/to/prodigal_genes.fna"
    metadata_file: "/path/to/metadata.txt"
    gtdb_file: "/path/to/gtdb_taxonomy.tsv"
    
    # Outputs
    root_out: "/path/to/output_directory"
    
    # Parameters
    timepoints_combinations:
      - ["pre", "post"]
    
    groups_combinations:
      - ["control", "treatment"]
    
    # Analysis options
    analysis_options:
      use_lmm: True
      use_significance_tests: True
      data_type: "longitudinal"  # or "single"
    
    min_sample_num: 4
    breadth_threshold: 0.1
    disable_zero_diff_filtering: False
    alpha: 0.05
    test_type: "both"
    preprocess_two_sample: True