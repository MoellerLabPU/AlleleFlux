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
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        allele_freq=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_single.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_mean.tsv.gz"
            )
        ),

        allele_freq_no_zero_diff=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_no_constant.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_no_zero-diff.tsv.gz"
            )
        ) if (
            (DATA_TYPE == "longitudinal" and not config["quality_control"].get("disable_zero_diff_filtering", False))
            or (DATA_TYPE == "single" and not config["quality_control"].get("disable_zero_diff_filtering", False))
        ) else [],

        longitudinal_output=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_longitudinal.tsv.gz"
        ) if DATA_TYPE == "longitudinal" else [],
    params:
        outDir=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
        ),
        fasta=config["input"]["fasta_path"],
        breath_threshold=config["quality_control"].get("breadth_threshold", 0.1),
        disable_zero_diff_filtering=(
            "--disable_zero_diff_filtering"
            if config["quality_control"].get("disable_zero_diff_filtering", False)
            else ""
        ),
        # Use the global DATA_TYPE variable
        data_type=DATA_TYPE,        

    threads: config["resources"]["cpus"]["analyze_alleles"]
    resources:
        mem_mb=config["resources"]["memory"]["analyze_alleles"],
        time=config["resources"]["time"]["general"],
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
            {params.disable_zero_diff_filtering} \
            --mag_mapping_file {input.mag_mapping}
            """


