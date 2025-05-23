Tutorial
========

This tutorial demonstrates how to analyze a simple dataset with AlleleFlux.

Prerequisites
---------------

Before starting this tutorial, make sure you have:

1. Installed AlleleFlux (see :doc:`../getting_started/installation`)
2. Downloaded the example data (if available)
3. Basic familiarity with the command line

Step 1: Prepare Input Files
----------------------------

First, prepare your input files:

1. **BAM Files**: Aligned reads for each sample
2. **Reference Genome**: FASTA file with reference sequences
3. **Prodigal Genes**: Predicted genes in the reference genome
4. **Metadata File**: Sample metadata with BAM paths

Example metadata file (``metadata.tsv``):

.. code-block:: text

    sample_id    bam_path                               subjectID    group      replicate    time
    S1           /path/to/bamfiles/S1.sorted.bam        mouse1       control    A           pre
    S2           /path/to/bamfiles/S2.sorted.bam        mouse2       control    B           pre
    S3           /path/to/bamfiles/S3.sorted.bam        mouse3       treatment  A           pre
    S4           /path/to/bamfiles/S4.sorted.bam        mouse4       treatment  B           pre
    S5           /path/to/bamfiles/S5.sorted.bam        mouse1       control    A           post
    S6           /path/to/bamfiles/S6.sorted.bam        mouse2       control    B           post
    S7           /path/to/bamfiles/S7.sorted.bam        mouse3       treatment  A           post
    S8           /path/to/bamfiles/S8.sorted.bam        mouse4       treatment  B           post

**Note**: The metadata file must include a ``bam_path`` column with full paths to BAM files. You can use the ``alleleflux-add-bam-path`` utility to add this column automatically.

Step 2: Configure the Pipeline
-------------------------------

Create a configuration file (``config.yml``) based on the template:

.. code-block:: yaml

    # Data type
    data_type: "longitudinal"  # or "single" for single timepoint

    # Input files
    input:
      bam_dir: "/path/to/bamfiles"  # For backward compatibility
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
      breadth_threshold: 0.1
      disable_zero_diff_filtering: false
    
    # Analysis options
    analysis:
      use_lmm: true
      use_significance_tests: true
      use_cmh: true  # Enable Cochran-Mantel-Haenszel tests
      significance_threshold: 0.05

Step 3: Run the Pipeline
------------------------

Method 1: Using the Main Pipeline Runner (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Run the complete workflow
    python alleleFlux.py --config config.yml

    # Or run with specific options
    python alleleFlux.py --config config.yml --threads 16 --profile cornell_profile/

Method 2: Using Snakemake Directly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Step 1: Profile samples and generate eligibility table
    cd smk_workflow
    snakemake -s step1.smk --configfile ../config.yml
    
    # Step 2: Analyze alleles and calculate scores
    snakemake -s step2.smk --configfile ../config.yml

Step 4: Examine the Results
----------------------------

1. **Eligibility Table**:
   
   Check which MAGs have sufficient coverage:
   
   .. code-block:: bash
   
       cat /path/to/output_directory/eligibility_table_pre_post-control_treatment.tsv

2. **Allele Frequency Analysis**:
   
   Check the allele frequency analysis results:
   
   .. code-block:: bash
   
       ls /path/to/output_directory/allele_analysis/allele_analysis_pre_post-control_treatment/

3. **Statistical Test Results**:
   
   Look at the statistical test results:
   
   .. code-block:: bash
   
       ls /path/to/output_directory/significance_tests/

4. **Scores**:
   
   Look at the MAG scores:
   
   .. code-block:: bash
   
       cat /path/to/output_directory/scores/processed/combined/scores_two_sample_unpaired-pre_post-control_treatment-MAGs.tsv

5. **Outlier Genes**:
   
   Identify genes with strong selection:
   
   .. code-block:: bash
   
       ls /path/to/output_directory/outlier_genes/pre_post-control_treatment/

Step 5: Understanding the Output
-------------------------------

The AlleleFlux workflow generates several types of output:

**Parallelism Scores**: Measure how consistently allele frequencies change across replicates within each group.

**Divergence Scores**: Quantify the degree of allele-frequency divergence between experimental groups.

**Statistical Tests**: Results from various statistical approaches (two-sample, single-sample, LMM, CMH).

**Outlier Genes**: Genes with exceptionally high scores that may be under strong selection.

Step 6: Advanced Analysis
------------------------

For more detailed analysis, you can use individual command-line tools:

.. code-block:: bash

    # Run CMH test separately
    alleleflux-cmh --input_df /path/to/longitudinal.tsv.gz \
        --preprocessed_df /path/to/preprocessed.tsv.gz \
        --min_sample_num 4 --mag_id MAG_ID --data_type longitudinal \
        --cpus 16 --output_dir /path/to/output

    # Calculate gene scores
    alleleflux-gene-scores --scores_file /path/to/scores.tsv \
        --output_dir /path/to/gene_scores

    # Detect outlier genes
    alleleflux-outliers --scores_file /path/to/gene_scores.tsv \
        --output_dir /path/to/outliers

Step 7: Visualize the Results
-----------------------------

You can visualize the results using your favorite plotting tools (e.g., R, Python). For example, to create a simple plot of scores across MAGs:

.. code-block:: python

    import pandas as pd
    import matplotlib.pyplot as plt
    
    # Load scores
    scores = pd.read_csv("scores_two_sample_unpaired-pre_post-control_treatment-MAGs.tsv", sep="\t")
    
    # Plot parallelism scores
    plt.figure(figsize=(10, 6))
    plt.bar(scores["MAG_ID"], scores["score_two_sample_unpaired (%)"])
    plt.xlabel("MAG ID")
    plt.ylabel("Parallelism Score (%)")
    plt.title("Parallelism Scores Across MAGs")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

Next Steps
----------

For more advanced usage and detailed explanations, see:

* :doc:`../usage/running_workflow` - Detailed workflow documentation
* :doc:`../usage/interpreting_results` - How to interpret AlleleFlux outputs
* :doc:`use_cases` - Real-world application examples
    plt.savefig("parallelism_scores.png")

Conclusion
-----------

In this tutorial, we've demonstrated how to:

1. Prepare input files for AlleleFlux
2. Configure the pipeline
3. Run the Snakemake workflow
4. Examine the results
5. Create simple visualizations

Next Steps
-----------

- Try adjusting parameters in the configuration file to see how they affect the results
- Apply this workflow to your own data
- Explore more advanced visualizations and statistical analyses