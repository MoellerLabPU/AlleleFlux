"""
Workflow for Step 2: Allele Analysis and Scoring
This workflow is designed to analyze alleles and score them based on various criteria.
"""

# Include modular workflow components
include: "rules/common.smk"
include: "rules/targets.smk"
include: "rules/allele_analysis.smk"
include: "rules/significance_tests.smk"
include: "rules/scoring.smk"
include: "rules/gene_analysis.smk"

# Define the main rule for this workflow
localrules:
    all,

rule all:
    input:
        get_allele_analysis_targets(),
        # lambda wildcards: get_two_sample_targets("two_sample_unpaired") if config["analysis_options"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_two_sample_targets("two_sample_paired") if config["analysis_options"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_single_sample_targets() if config["analysis_options"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_lmm_targets() if config["analysis_options"].get("use_lmm", True) else [],
        # get_significance_scores_targets(),
        # get_combined_scores_targets(),
        get_taxa_scores_targets(),
        # get_gene_scores_targets(),
        get_outlier_gene_targets(),
