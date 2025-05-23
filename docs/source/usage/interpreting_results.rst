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
   - ``{mag}_cmh.tsv.gz``: Results of Cochran-Mantel-Haenszel tests

   Each file contains p-values and test statistics for genomic positions. Key columns typically include:
   
   - **contig**: Contig identifier
   - **position**: Genomic position (0-based)
   - **gene_id**: Gene identifier (when available)
   - **p_value_{test_type}**: P-value from the statistical test
   - Additional test-specific columns (e.g., effect sizes, confidence intervals)

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

Understanding Output File Formats
----------------------------------

AlleleFlux generates several types of output files with specific formats:

1. **Profile Files**

   Found in ``profiles/{sample}.sorted/{sample}_{mag}_profiled.tsv.gz``:
   
   - **contig**: Contig identifier
   - **position**: Genomic position (0-based)
   - **ref_base**: Reference base at this position
   - **total_coverage**: Total read coverage
   - **A, C, G, T, N**: Counts for each nucleotide
   - **gene_id**: Gene identifier (when mapped to genes)

2. **Allele Frequency Files**

   Found in ``allele_analysis/allele_analysis_{timepoints}-{groups}/{mag}_allele_frequency_*.tsv.gz``:
   
   - Contains processed allele frequency data
   - Includes sample metadata (group, timepoint, replicate)
   - May be filtered based on quality control settings

3. **Gene Score Files**

   Gene-level score files contain detailed information about individual genes:
   
   - **gene_id**: Gene identifier
   - **total_sites_per_group_{test}**: Total sites tested for this gene
   - **significant_sites_per_group_{test}**: Number of significant sites
   - **score_{test} (%)**: Percentage score for this gene

4. **Outlier Gene Files**

   Outlier detection results include statistical assessments:
   
   - **gene_id**: Gene identifier
   - **mag_score_{test} (%)**: MAG-level score for comparison
   - **gene_score_{test} (%)**: Gene-level score
   - **p_value_binomial_{test}**: Binomial test p-value
   - **p_value_poisson_{test}**: Poisson test p-value

Understanding CMH Test Results
-----------------------------

The Cochran-Mantel-Haenszel (CMH) test provides a powerful approach to identifying significant allele frequency changes across multiple timepoints while controlling for confounding factors.

1. **CMH Test Output Files**

   The main output files for CMH tests are found in:
   
   .. code-block:: text
   
       significance_tests/cmh_{timepoints}-{groups}/{mag}_cmh.tsv.gz
   
   These files contain:
   
   - **mag_id**: MAG identifier
   - **contig**: Contig identifier
   - **gene_id**: Gene identifier (if applicable)
   - **position**: Position on the contig (0-based)
   - **num_pairs**: Number of replicate pairs used in the test
   - **p_value_CMH**: P-value from the CMH test
   - **time**: Time point identifier (for longitudinal data)
   - **notes**: Any error messages or warnings

2. **CMH Scores**

   CMH scores are calculated based on the significance of allele frequency changes between timepoints. These scores are found in:
   
   .. code-block:: text
   
       scores/intermediate/MAG_scores_cmh_{timepoints}-{groups}/{mag}_score_cmh_{focus}.tsv
   
   Key columns include:
   
   - **MAG_ID**: MAG identifier
   - **focus_timepoint**: The focus timepoint used for analysis
   - **significant_sites_per_group_CMH**: Number of sites with significant changes
   - **total_sites_per_group_CMH**: Total number of tested sites
   - **score_CMH (%)**: Percentage of sites with significant changes
   - **grouped_by**: Grouping method used (typically "MAG_ID")

3. **Interpreting CMH Results**

   - **High CMH scores** indicate consistent allele frequency changes across replicates between timepoints
   - **Focus timepoint significance** helps identify timepoint-specific selection events  
   - **Comparison across MAGs/genes** reveals differences in selective pressure
   - **CMH tests support three analysis modes**: single timepoint, longitudinal, and across-time comparisons

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

Here's an example workflow for interpreting AlleleFlux results:

1. **Start with Eligibility Analysis**
   
   Examine the eligibility table to see which MAGs have sufficient coverage:
   
   .. code-block:: bash
   
       head eligibility_table_pre_post-control_treatment.tsv

2. **Examine MAG-level Scores**
   
   Look at the combined MAG scores to identify MAGs with high parallelism or divergence:
   
   .. code-block:: bash
   
       # Two-sample unpaired test results
       head scores_two_sample_unpaired-pre_post-control_treatment-MAGs.tsv
       
       # CMH test results (if available)
       head scores_cmh-pre_post-control_treatment-MAGs.tsv

3. **Analyze Taxonomic Patterns**
   
   Examine taxonomic aggregations to identify groups showing strong signals:
   
   .. code-block:: bash
   
       # Family-level aggregation
       head scores_two_sample_unpaired-pre_post-control_treatment-family.tsv

4. **Investigate Gene-level Results**
   
   Focus on specific MAGs with high scores and examine gene-level results:
   
   .. code-block:: bash
   
       # Individual gene scores for a specific MAG
       head MAG123_two_sample_unpaired_gene_scores_individual.tsv
       
       # Outlier genes for the same MAG
       head MAG123_two_sample_unpaired_outlier_genes.tsv

5. **Compare Statistical Approaches**
   
   Compare results from different statistical tests for consistency:
   
   - Two-sample tests (paired vs unpaired)
   - Linear Mixed Models (LMM)
   - Cochran-Mantel-Haenszel (CMH) tests

6. **Validate Interesting Findings**
   
   For genes identified as outliers or showing high scores:
   
   - Check their functional annotations
   - Examine their genomic context
   - Consider their biological relevance to your experimental conditions

Further Analysis
-----------------

After identifying genes of interest, you may want to:

1. Investigate the specific mutations within these genes
2. Compare the results across different statistical tests
3. Correlate the findings with metadata (e.g., clinical outcomes, environmental factors)
4. Validate key findings with follow-up experiments or targeted sequencing

Troubleshooting Interpretation Issues
-------------------------------------

**Empty or Missing Results**

If you encounter empty result files:

- Check the eligibility table to ensure MAGs meet minimum sample requirements
- Verify that input data has sufficient coverage and quality
- Review quality control parameters (``breadth_threshold``, ``min_sample_num``)

**Low Scores Across All MAGs**

If all scores are consistently low:

- Consider adjusting the p-value threshold in configuration
- Check if experimental conditions provide sufficient selective pressure
- Verify that timepoints are appropriate for detecting evolutionary changes

**Inconsistent Results Between Statistical Tests**

If different tests give conflicting results:

- LMM tests may be more sensitive to experimental design complexities
- Two-sample tests may be affected by unbalanced group sizes
- CMH tests are particularly suited for detecting consistent directional changes

**High Background Significance**

If many positions show significance without biological relevance:

- Consider more stringent p-value thresholds
- Enable preprocessing to filter out low-effect sites
- Check for technical artifacts in sequencing or alignment

**Missing Gene Annotations**

If gene IDs are missing from results:

- Ensure Prodigal gene prediction was run on the reference sequences
- Verify that the gene FASTA file path is correct in configuration
- Check that contig names match between reference FASTA and gene predictions