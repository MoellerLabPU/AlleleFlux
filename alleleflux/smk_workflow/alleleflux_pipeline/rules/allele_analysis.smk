"""Allele frequency analysis rules.

This module contains rules for analyzing allele frequencies from MAG profiles.
It generates allele frequency tables and optionally filters positions with zero
frequency differences (constant positions) across samples.
"""

rule analyze_alleles:
    input:
        qc_file=os.path.join(
            OUTDIR,
            "QC",
            "QC_{timepoints}-{groups}",
            "{mag}_QC.tsv",
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
                "{mag}_allele_frequency_single.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_mean.tsv.gz"
            )
        ),
        # Simplified conditional: disable_zero_diff_filtering check is sufficient
        # since DATA_TYPE is already evaluated at module load time
        allele_freq_no_zero_diff=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_no_constant.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_no_zero-diff.tsv.gz"
            )
        ) if not config["quality_control"].get("disable_zero_diff_filtering", False) else [],

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
        disable_zero_diff_filtering=(
            "--disable_zero_diff_filtering"
            if config["quality_control"].get("disable_zero_diff_filtering", False)
            else ""
        ),
        # Use the global DATA_TYPE variable
        data_type=DATA_TYPE,

    threads: get_threads("allele_analysis")
    resources:
        mem_mb=get_mem_mb("allele_analysis"),
        time=get_time("allele_analysis"),
    shell:
        """
        alleleflux-allele-freq \
            --magID {wildcards.mag} \
            --qc_file {input.qc_file} \
            --data_type {params.data_type} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.disable_zero_diff_filtering}
        """


