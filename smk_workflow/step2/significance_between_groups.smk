def get_between_group_inputs(test_type=None):
    if test_type == "cmh" and DATA_TYPE == "longitudinal":
        # For CMH test, use the preprocessed file
        between_groups_input = os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_longitudinal.tsv.gz"
        )
        return between_groups_input

    preprocess_enabled = config["statistics"].get("preprocess_two_sample", False)
    if preprocess_enabled:
        if DATA_TYPE == "single":
            between_groups_input = os.path.join(
                OUTDIR,
                "significance_tests",
                "preprocessed_between_groups_{timepoints}-{groups}",
                "{mag}_allele_frequency_preprocessed.tsv.gz",
            )
        elif DATA_TYPE == "longitudinal":
            between_groups_input = os.path.join(
                OUTDIR,
                "significance_tests",
                "preprocessed_between_groups_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz",
            )
    else:
        # Input file depends on data_type and filtering options
        if DATA_TYPE == "single":
            if not config["quality_control"].get("disable_zero_diff_filtering", False):
                # When single data type and filtering is not disabled
                between_groups_input = os.path.join(
                    OUTDIR,
                    "allele_analysis",
                    "allele_analysis_{timepoints}-{groups}",
                    "{mag}_allele_frequency_no_constant.tsv.gz",
                )
            else:
                # When single data type and filtering is disabled
                between_groups_input = os.path.join(
                    OUTDIR,
                    "allele_analysis",
                    "allele_analysis_{timepoints}-{groups}",
                    "{mag}_allele_frequency_single.tsv.gz",
                )
        elif DATA_TYPE == "longitudinal":
            between_groups_input = os.path.join(
                OUTDIR,
                "allele_analysis",
                "allele_analysis_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean.tsv.gz",
            )
    return between_groups_input

rule preprocess_between_groups:
    input:
        os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_no_constant.tsv.gz"
                if (
                    DATA_TYPE == "single"
                    and not config["quality_control"].get("disable_zero_diff_filtering", False)
                )
                else (
                    "{mag}_allele_frequency_single.tsv.gz"
                    if DATA_TYPE == "single"
                    else "{mag}_allele_frequency_changes_mean.tsv.gz"
                )
            ),
        ),
    output:
        outPath=os.path.join(
            OUTDIR,
            "significance_tests",
            "preprocessed_between_groups_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_preprocessed.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz"
            ),
        ),
    params:
        alpha=config["statistics"].get("alpha", 0.05),
        filter_type=config["statistics"].get("filter_type", "t-test"),
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["quality_control"]
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-preprocess-between-groups \
            --mean_changes_fPath {input} \
            --cpus {threads} --alpha {params.alpha} \
            --output_fPath {output.outPath} \
            --filter_type {params.filter_type} \
            --data_type {params.data_type}
            """
            
rule two_sample_unpaired:
    input:
        input_df=get_between_group_inputs(),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_unpaired_{timepoints}-{groups}",
            "{mag}_two_sample_unpaired.tsv.gz",
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_unpaired_{timepoints}-{groups}"
        ),
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["significance_test"]
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-two-sample-unpaired \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type}
            """

rule two_sample_paired:
    input:
        input_df=get_between_group_inputs(),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_paired_{timepoints}-{groups}",
            "{mag}_two_sample_paired.tsv.gz",
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_paired_{timepoints}-{groups}"
        ),
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["significance_test"]
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-two-sample-paired \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type}
            """

rule lmm_analysis:
    input:
        input_df=get_between_group_inputs(),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "lmm_{timepoints}-{groups}",
            "{mag}_lmm.tsv.gz",
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "lmm_{timepoints}-{groups}"
        ),
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["significance_test"]
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-lmm \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type} \
            --mag_id {wildcards.mag}
            """

rule cmh_test:
    input:
        input_df=get_between_group_inputs(test_type="cmh"),
        preprocessed_df=os.path.join(
                OUTDIR,
                "significance_tests",
                "preprocessed_between_groups_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz"
            )
            if config["statistics"].get("preprocess_two_sample", True) and DATA_TYPE == "longitudinal"
            else [],
    
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "cmh_{timepoints}-{groups}",
            "{mag}_cmh.tsv.gz"
        )
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "cmh_{timepoints}-{groups}"
        ),
        data_type=DATA_TYPE,
        # Conditionally include the preprocessed file argument.
        preprocessed_flag=(
            ("--preprocessed_df")
            if config["statistics"].get("preprocess_two_sample", True) and DATA_TYPE == "longitudinal"
            else ""
        ),
    threads: config["resources"]["cpus"]["significance_test"]
    resources:
        time=config["resources"]["time"]["significance_test"],
        mem_mb=config["resources"]["memory"]["significance_test"],
    shell:
        """
        alleleflux-cmh \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --data_type {params.data_type} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.preprocessed_flag} {input.preprocessed_df}
            
        """

