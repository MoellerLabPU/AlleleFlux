def get_two_sample_inputs():
    preprocess_enabled = config["statistics"].get("preprocess_two_sample", False)
    if preprocess_enabled:
        if DATA_TYPE == "single":
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
        if DATA_TYPE == "single":
            if not config["quality_control"].get("disable_zero_diff_filtering", False):
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

rule preprocess_two_sample:
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
            "preprocessed_two_sample_{timepoints}-{groups}",
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
        alleleflux-preprocess-two-sample \
            --mean_changes_fPath {input} \
            --cpus {threads} --alpha {params.alpha} \
            --output_fPath {output.outPath} \
            --filter_type {params.filter_type} \
            --data_type {params.data_type}
            """
            
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
        input_df=get_two_sample_inputs(),
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
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "single_sample_{timepoints}-{groups}"
        ),
        max_zero_flag=(
            "--max_zero_count " + str(config["statistics"]["max_zero_count"])
            if config["statistics"].get("max_zero_count", None) is not None
            else ""
        ),
    threads: config["resources"]["cpus"]["significance_test"]
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-single-sample \
            --mean_changes_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --group {wildcards.group} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.max_zero_flag}
            """

rule lmm_analysis:
    input:
        input_df=get_two_sample_inputs(),
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
