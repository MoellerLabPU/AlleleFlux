def get_two_sample_inputs():
    preprocess_enabled = config.get("preprocess_two_sample", False)
    if preprocess_enabled:
        if config["data_type"] == "single":
            two_sample_input = os.path.join(
                OUTDIR,
                "significance_tests",
                "preprocessed_two_sample_{timepoints}-{groups}",
                "{mag}_allele_frequency_preprocessed.tsv.gz",
            )
        else:  # longitudinal
            two_sample_input = os.path.join(
                OUTDIR,
                "significance_tests",
                "preprocessed_two_sample_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz",
            )
    else:
        # Input file depends on data_type and filtering options
        if config["data_type"] == "single":
            if not config.get("disable_zero_diff_filtering", False):
                # When single data type and filtering is not disabled
                two_sample_input = os.path.join(
                    OUTDIR,
                    "allele_analysis",
                    "allele_analysis_{timepoints}-{groups}",
                    "{mag}_allele_frequency_no_constant.tsv.gz",
                )
            else:
                # When single data type and filtering is disabled
                two_sample_input = os.path.join(
                    OUTDIR,
                    "allele_analysis",
                    "allele_analysis_{timepoints}-{groups}",
                    "{mag}_allele_frequency_single.tsv.gz",
                )
        else:  # longitudinal
            two_sample_input = os.path.join(
                OUTDIR,
                "allele_analysis",
                "allele_analysis_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean.tsv.gz",
            )
    return two_sample_input


rule two_sample_unpaired:
    input:
        input_df=get_two_sample_inputs(),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_unpaired_{timepoints}-{groups}",
            "{mag}_two_sample_unpaired.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_unpaired_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["two_sample_unpaired"],
        data_type=config["data_type"],
    threads: config["cpus"]["significance_test"]
    resources:
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type}
            """


rule two_sample_paired:
    input:
        input_df=get_two_sample_inputs(),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_paired_{timepoints}-{groups}",
            "{mag}_two_sample_paired.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_paired_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["two_sample_paired"],
        data_type=config["data_type"],
    threads: config["cpus"]["significance_test"]
    resources:
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --input_df {input.input_df} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type}
            """


rule single_sample:
    input:
        mean_allele_changes=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean.tsv.gz",
        ),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "single_sample_{timepoints}-{groups}",
            "{mag}_single_sample_{group}.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "single_sample_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["single_sample"],
        max_zero_flag=(
            "--max_zero_count " + str(config["max_zero_count"])
            if config.get("max_zero_count", None) is not None
            else ""
        ),
    threads: config["cpus"]["significance_test"]
    resources:
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --mean_changes_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --group {wildcards.group} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.max_zero_flag}
            """
