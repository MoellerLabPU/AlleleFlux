"""
Workflow for Step 2: Allele Analysis and Scoring
This workflow is designed to analyze alleles and score them based on various criteria.
"""


# Include modular workflow components
include: "common.smk"
include: "targets.smk"
include: "allele_analysis.smk"
include: "significance_tests.smk"
include: "scoring.smk"
include: "gene_analysis.smk"


# Define the main rule for this workflow
localrules:
    all,


rule all:
    input:
        get_allele_analysis_targets(),
        # get_two_sample_targets("two_sample_unpaired"),
        # get_two_sample_targets("two_sample_paired"),
        # get_single_sample_targets(),
        # get_significance_scores_targets(),
        # get_combined_scores_targets(),
        get_taxa_scores_targets(),
        # get_gene_scores_targets(),
        get_outlier_gene_targets(),
