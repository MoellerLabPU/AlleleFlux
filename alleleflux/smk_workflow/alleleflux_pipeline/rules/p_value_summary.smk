"""P-value summary aggregation rules.

This module aggregates p-values from significance tests across all MAGs and
positions, applying FDR correction (if configured) and generating summary
tables for downstream analysis including dN/dS calculations.
"""


def get_all_significance_test_results_for_summary(wildcards):
    """
    Get significance results for a specific test type, used as input for p_value_summary.
    
    Respects preprocessing eligibility by triggering the appropriate checkpoint
    and using get_eligible_mags() instead of direct eligibility lookup.
    """
    targets = []
    test_type = wildcards.test_type
    timepoints = wildcards.timepoints
    groups = wildcards.groups
    
    # Check preprocessing config
    preprocess_between = config["statistics"].get("preprocess_between_groups", False)
    preprocess_within = config["statistics"].get("preprocess_within_groups", False)

    if test_type in ['two_sample_unpaired', 'two_sample_paired', 'lmm', 'cmh']:
        # Trigger preprocessing eligibility checkpoint if enabled
        if preprocess_between:
            checkpoints.preprocessing_eligibility_between_groups.get(
                timepoints=timepoints,
                groups=groups
            )
        
        mags = get_eligible_mags(timepoints, groups, test_type)
        for mag in mags:
            targets.append(
                os.path.join(
                    OUTDIR,
                    "significance_tests",
                    f"{test_type}_{timepoints}-{groups}",
                    f"{mag}_{test_type}.tsv.gz",
                )
            )

    elif test_type in ["single_sample",'lmm_across_time', 'cmh_across_time'] and DATA_TYPE == "longitudinal":
        # Trigger preprocessing eligibility checkpoint if enabled
        if preprocess_within:
            checkpoints.preprocessing_eligibility_within_groups.get(
                timepoints=timepoints,
                groups=groups
            )
        
        sample_entries = get_eligible_mags(timepoints, groups, test_type)
        for mag, group in sample_entries:
            targets.append(
                os.path.join(
                    OUTDIR,
                    "significance_tests",
                    f"{test_type}_{timepoints}-{groups}",
                    f"{mag}_{test_type}_{group}.tsv.gz",
                )
            )
    
    return targets


rule p_value_summary:
    input:
        get_all_significance_test_results_for_summary
    output:
        os.path.join(
            OUTDIR, 
            "p_value_summary", 
            "{timepoints}-{groups}", 
            "p_value_summary_{test_type}_{timepoints}.tsv"
        )
    params:
        output_dir=os.path.join(OUTDIR, "p_value_summary", "{timepoints}-{groups}"),
        input_dir=os.path.join(OUTDIR, "significance_tests"),
        prefix="p_value_summary",
        group_by_mag_id=(
            "--fdr-group-by-mag-id"
            if config["statistics"].get("fdr_group_by_mag_id", False)
            else ""
        )
    shell:
        """
        alleleflux-p-value-summary \
            --input-dir {params.input_dir} \
            --timepoints {wildcards.timepoints} \
            --test-types {wildcards.test_type} \
            --outdir {params.output_dir} \
            --prefix {params.prefix} \
            {params.group_by_mag_id}
        """
