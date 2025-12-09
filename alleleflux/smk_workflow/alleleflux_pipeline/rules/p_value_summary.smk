def get_all_significance_test_results_for_summary(wildcards):
    """
    Get significance results for a specific test type, used as input for p_value_summary.
    """
    targets = []
    test_type = wildcards.test_type
    timepoints = wildcards.timepoints
    groups = wildcards.groups

    if test_type in ['two_sample_unpaired', 'two_sample_paired', 'lmm', 'cmh']:
        mags = get_mags_by_eligibility(
            timepoints, groups, eligibility_type=test_type
        )
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
        sample_entries = get_single_sample_entries(timepoints, groups)
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
