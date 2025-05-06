Interpreting Results
=====================

This guide explains how to interpret the results produced by AlleleFlux.

Output Directory Structure
---------------------------

After running the workflow, the output directory will have the following structure:

.. code-block:: text

    root_out/
    ├── single_timepoint/ (or longitudinal/)
        ├── profiles/
        │   └── {sample}.sorted/
        │       └── {sample}_{mag}_profiled.tsv.gz
        ├── inputMetadata/
        │   └── inputMetadata_{timepoints}-{groups}/
        │       └── {mag}_metadata.tsv
        ├── QC/
        │   └── QC_{timepoints}-{groups}/
        │       └── {mag}_qc.tsv
        ├── eligibility_table_{timepoints}-{groups}.tsv
        ├── allele_analysis/
        │   └── allele_analysis_{timepoints}-{groups}/
        │       └── {mag}_allele_frequency_*.tsv.gz
        ├── significance_tests/
        │   ├── two_sample_unpaired_{timepoints}-{groups}/
        │   │   └── {mag}_two_sample_unpaired.tsv.gz
        │   ├── two_sample_paired_{timepoints}-{groups}/
        │   │   └── {mag}_two_sample_paired.tsv.gz
        │   ├── single_sample_{timepoints}-{groups}/
        │   │   └── {mag}_single_sample_{group}.tsv.gz
        │   ├── lmm_{timepoints}-{groups}/
        │   │   └── {mag}_lmm.tsv.gz
        │   └── cmh_{timepoints}-{groups}/
        │       └── {mag}_cmh.tsv.gz
        ├── scores/
        │   ├── intermediate/
        │   │   ├── MAG_scores_{timepoints}-{groups}/
        │   │   │   └── {mag}_score_{test_type}.tsv
        │   │   └── MAG_scores_cmh_{timepoints}-{groups}/
        │   │       └── {mag}_score_cmh.tsv
        │   └── processed/
        │       ├── combined/
        │       │   ├── scores_{test_type}-{timepoints}-{groups}-MAGs.tsv
        │       │   └── scores_{test_type}-{timepoints}-{groups}-{taxon}.tsv
        │       └── gene_scores_{timepoints}-{groups}/
        │           ├── {mag}_{test_type}_gene_scores_combined.tsv
        │           ├── {mag}_{test_type}_gene_scores_individual.tsv
        │           └── {mag}_{test_type}_gene_scores_overlapping.tsv
        └── outlier_genes/
            └── {timepoints}-{groups}/
                └── {mag}_{test_type}_outlier_genes.tsv

Key Result Files
-----------------

1. **Eligibility Table**

   ``eligibility_table_{timepoints}-{groups}.tsv``

   This file lists all MAGs and indicates which statistical tests they are eligible for based on coverage and sample requirements.

2. **Statistical Test Results**

   - ``{mag}_two_sample_unpaired.tsv.gz``: Results of unpaired two-sample tests
   - ``{mag}_two_sample_paired.tsv.gz``: Results of paired two-sample tests
   - ``{mag}_single_sample_{group}.tsv.gz``: Results of single-sample tests
   - ``{mag}_lmm.tsv.gz``: Results of linear mixed model tests

3. **Score Files**

   - ``{mag}_score_{test_type}.tsv``: Scores for each MAG
   - ``scores_{test_type}-{timepoints}-{groups}-MAGs.tsv``: Combined scores for all MAGs
   - ``scores_{test_type}-{timepoints}-{groups}-{taxon}.tsv``: Scores aggregated by taxonomic rank

4. **Gene Score Files**

   - ``{mag}_{test_type}_gene_scores_combined.tsv``: Combined gene scores
   - ``{mag}_{test_type}_gene_scores_individual.tsv``: Individual gene scores
   - ``{mag}_{test_type}_gene_scores_overlapping.tsv``: Scores for overlapping genes

5. **Outlier Genes**

   ``{mag}_{test_type}_outlier_genes.tsv``: Genes identified as outliers

Understanding Scores
---------------------

AlleleFlux calculates two key scores:

1. **Parallelism Score**

   The parallelism score quantifies the degree to which allele frequencies change in parallel across replicates within a group. High scores indicate strong parallel evolution, suggesting deterministic rather than stochastic changes.

   - Score range: 0-100%
   - Higher scores indicate stronger parallelism (more deterministic changes)
   - Significance is determined by p-values from statistical tests

2. **Divergence Score**

   The divergence score quantifies the degree to which allele frequencies diverge between groups. High scores indicate strong differential selection between groups.

   - Score range: 0-100%
   - Higher scores indicate stronger divergence between groups
   - Significance is determined by p-values from statistical tests

Understanding CMH Test Results
-----------------------------

The Cochran-Mantel-Haenszel (CMH) test provides a powerful approach to identifying significant allele frequency changes across multiple timepoints while controlling for confounding factors.

1. **CMH Test Output Files**

   The main output files for CMH tests are found in:
   
   .. code-block:: text
   
       significance_tests/cmh_{timepoints}-{groups}/{mag}_cmh.tsv.gz
   
   These files contain:
   
   - **contig_id**: Contig identifier
   - **position**: Position on the contig
   - **pvalue**: P-value for the CMH test
   - **time**: Time point identifier (if applicable)

2. **CMH Scores**

   CMH scores are calculated based on the significance of allele frequency changes between timepoints. These scores are found in:
   
   .. code-block:: text
   
       scores/intermediate/MAG_scores_cmh_{timepoints}-{groups}/{mag}_score_cmh.tsv
   
   Key columns include:
   
   - **mag_id**: MAG identifier
   - **significant_sites**: Number of sites with significant changes
   - **total_sites**: Total number of tested sites
   - **cmh_score**: Ratio of significant sites to total sites
   - **focus_timepoint**: The focus timepoint used for analysis

3. **Interpreting CMH Results**

   - **High CMH scores** indicate consistent allele frequency changes across replicates between timepoints
   - **Focus timepoint significance** helps identify timepoint-specific selection events
   - **Comparison across MAGs/genes** reveals differences in selective pressure

Identifying Significant Results
--------------------------------

To identify significant results:

1. **MAG-level Analysis**
   
   - Look for MAGs with high scores in ``scores_{test_type}-{timepoints}-{groups}-MAGs.tsv``
   - Focus on MAGs with high percentages of significant sites
   - Compare scores across different test types (paired, unpaired, LMM, CMH)

2. **Taxonomic Analysis**
   
   - Compare scores across taxonomic ranks in ``scores_{test_type}-{timepoints}-{groups}-{taxon}.tsv``
   - Look for taxonomic groups with consistently high scores
   - Consider scores from multiple statistical approaches for robust interpretation

3. **Gene-level Analysis**
   
   - Identify genes with exceptionally high scores in ``{mag}_{test_type}_gene_scores_individual.tsv``
   - Focus on outlier genes identified in ``{mag}_{test_type}_outlier_genes.tsv``
   - Pay special attention to genes that appear as outliers in multiple tests

Example Interpretation
----------------------

Here's an example of how to interpret the results:

1. Start by examining the eligibility table to see which MAGs have sufficient coverage for analysis
2. Look at the MAG scores to identify MAGs with high parallelism or divergence scores
3. Examine taxonomic scores to determine which taxonomic groups show strong signals
4. Compare results from different statistical tests (Two-sample, LMM, CMH) for consistency
5. Focus on outlier genes to identify specific genes under strong selection
6. For each outlier gene, investigate its function and potential relevance to the experimental conditions

Further Analysis
-----------------

After identifying genes of interest, you may want to:

1. Investigate the specific mutations within these genes
2. Compare the results across different statistical tests
3. Correlate the findings with metadata (e.g., clinical outcomes, environmental factors)
4. Validate key findings with follow-up experiments or targeted sequencing