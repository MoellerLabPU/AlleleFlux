"""
Workflow for Step 2: Allele Analysis and Scoring
This workflow is designed to analyze alleles and score them based on various criteria.
"""

# from alleleflux.scripts.utilities.logging_config import setup_logging
# setup_logging()

# Include modular workflow components
include: "shared/common.smk"

include: "step2/targets.smk"
include: "step2/allele_analysis.smk"
include: "step2/significance_within_group.smk"
include: "step2/significance_between_groups.smk"
include: "step2/scoring.smk"
include: "step2/gene_analysis.smk"
include: "step2/p_value_summary.smk"
include: "step2/dnds_analysis.smk"

# Define the main rule for this workflow
localrules:
    all,

rule all:
    input:
        get_allele_analysis_targets(),
        # Conditional targets based on allele_analysis_only setting
        lambda wildcards: [] if config["analysis"].get("allele_analysis_only", False) else get_taxa_scores_targets(),
        lambda wildcards: [] if config["analysis"].get("allele_analysis_only", False) else get_outlier_gene_targets(),
        lambda wildcards: [] if config["analysis"].get("allele_analysis_only", False) else get_p_value_summary_targets(),
        lambda wildcards: [] if config["analysis"].get("allele_analysis_only", False) else get_dnds_analysis_targets()
        # lambda wildcards: get_two_sample_targets("two_sample_unpaired") if config["analysis"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_two_sample_targets("two_sample_paired") if config["analysis"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_single_sample_targets() if config["analysis"].get("use_significance_tests", True) else [],
        # lambda wildcards: get_lmm_targets() if config["analysis"].get("use_lmm", True) else [],
        # get_significance_scores_targets(),
        # get_combined_scores_targets(),
        # Uncomment when the targets are defined
        # get_taxa_scores_targets(),
        # get_outlier_gene_targets(),

