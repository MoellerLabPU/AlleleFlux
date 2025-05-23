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

    # Data type: "single" for a single timepoint or "longitudinal" for multiple timepoints
    data_type: "longitudinal"

    # Input files
    input:
      bam_dir: "/path/to/bamfiles"  # Used for backward compatibility
      fasta_path: "/path/to/reference.fa"
      prodigal_path: "/path/to/prodigal_genes.fna"
      metadata_path: "/path/to/metadata.tsv"  # Must include bam_path column
      gtdb_path: "/path/to/gtdb_taxonomy.tsv"  # Optional
    
    # Output directory
    output:
      root_dir: "/path/to/output_directory"
    
    # Experimental design
    timepoints_combinations:
      - timepoint: ["pre", "post"]
        focus: "post"  # Required for CMH test with longitudinal data
    
    groups_combinations:
      - ["control", "treatment"]
    
    # Quality control parameters
    quality_control:
      min_sample_num: 4
      breadth_threshold: 0.1  # Note: 'breadth' not 'breath'
      disable_zero_diff_filtering: false
    
    # Analysis options
    analysis:
      use_lmm: true
      use_significance_tests: true
      use_cmh: true                # Enable Cochran-Mantel-Haenszel tests
      significance_threshold: 0.05  # Alpha value for statistical tests
    
    # Focus timepoints for CMH test (specify which timepoint to focus on)
    focus_timepoints:
      pre_post: "post"  # For CMH analysis focusing on the 'post' timepoint

Adding BAM Paths to Metadata
---------------------------

AlleleFlux includes a utility script to help you add BAM file paths to your existing metadata file:

.. code-block:: bash

    alleleflux-add-bam-path \
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

Creating MAG Mapping Files
--------------------------

Before running AlleleFlux, you may need to create a combined FASTA file and MAG mapping from individual MAG FASTA files:

.. code-block:: bash

    alleleflux-create-mag-mapping \
        --dir /path/to/mag_fastas \
        --extension fa \
        --output-fasta combined_mags.fasta \
        --output-mapping mag_mapping.tsv \
        --output-genomes-dir modified_genomes/

This preprocessing utility:

1. Combines individual MAG FASTA files into a single file
2. Creates a mapping file linking contigs to MAG IDs  
3. Optionally creates individual genome files with modified headers
4. Should be run BEFORE starting the main AlleleFlux workflow

Focus Timepoints for CMH Test
----------------------------

The Cochran-Mantel-Haenszel (CMH) test requires a focus timepoint to be specified for each timepoint combination. The focus timepoint is used to identify significant allele frequency changes between groups at a specific timepoint.

In the configuration file, specify the focus timepoint for each timepoint combination using the ``focus_timepoints`` section. The CMH test will calculate significance scores based on the specified focus timepoint, which is particularly useful for identifying evolutionary changes at specific points in time.