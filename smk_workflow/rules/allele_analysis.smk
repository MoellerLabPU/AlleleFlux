rule analyze_alleles:
    input:
        mag_metadata_file=os.path.join(
            OUTDIR,
            "inputMetadata",
            "inputMetadata_{timepoints}-{groups}",
            "{mag}_metadata.tsv",
        ),
        eligibility_table=os.path.join(
            OUTDIR,
            "eligibility_table_{timepoints}-{groups}.tsv",
        ),
    output:
        allele_freq=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_no_constant.tsv.gz"
                if (
                    DATA_TYPE == "single"
                    and not config.get("disable_zero_diff_filtering", False)
                )
                else (
                    "{mag}_allele_frequency_single.tsv.gz"
                    if DATA_TYPE == "single"
                    else "{mag}_allele_frequency_changes_mean.tsv.gz"
                )
            ),
        ),
    params:
        outDir=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
        ),
        fasta=config["fasta"],
        breath_threshold=config.get("breath_threshold", 0.1),
        disable_zero_diff_filtering=(
            "--disable_zero_diff_filtering"
            if config.get("disable_zero_diff_filtering", False)
            else ""
        ),
        # Use the global DATA_TYPE variable
        data_type=DATA_TYPE,
    threads: config["cpus"]["analyze_alleles"]
    resources:
        mem_mb=config["memory"]["analyze_alleles"],
        time=config["time"]["general"],
    shell:
        """
        alleleflux-allele-freq \
            --magID {wildcards.mag} \
            --mag_metadata_file {input.mag_metadata_file} \
            --fasta {params.fasta} \
            --breath_threshold {params.breath_threshold} \
            --data_type {params.data_type} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.disable_zero_diff_filtering}
            """

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
                    and not config.get("disable_zero_diff_filtering", False)
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
        alpha=config.get("alpha", 0.05),
        test_type=config.get("test_type", "t-test"),
        data_type=DATA_TYPE,
    threads: config["cpus"]["quality_control"]
    resources:
        time=config["time"]["general"],
    shell:
        """
        alleleflux-preprocess-two-sample \
            --mean_changes_fPath {input} \
            --cpus {threads} --alpha {params.alpha} \
            --output_fPath {output.outPath} \
            --test_type {params.test_type} \
            --data_type {params.data_type}
            """
