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
4. **Metadata File**: Sample metadata

Example metadata file (``metadata.txt``):

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

Step 2: Configure the Pipeline
-------------------------------

Create a configuration file (``config.yml``) based on the template:

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
      data_type: "longitudinal"
    
    min_sample_num: 4
    breadth_threshold: 0.1
    disable_zero_diff_filtering: False
    alpha: 0.05
    test_type: "both"
    preprocess_two_sample: True

Step 3: Run the Pipeline
------------------------

Run the Snakemake workflow:

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
   
       cat /path/to/output_directory/longitudinal/eligibility_table_pre_post-control_treatment.tsv

2. **Scores**:
   
   Look at the MAG scores:
   
   .. code-block:: bash
   
       cat /path/to/output_directory/longitudinal/scores/processed/combined/scores_two_sample_unpaired-pre_post-control_treatment-MAGs.tsv

3. **Outlier Genes**:
   
   Identify genes with strong selection:
   
   .. code-block:: bash
   
       cat /path/to/output_directory/longitudinal/outlier_genes/pre_post-control_treatment/MAG_ID_two_sample_unpaired_outlier_genes.tsv

Step 5: Visualize the Results
-------------------------------

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