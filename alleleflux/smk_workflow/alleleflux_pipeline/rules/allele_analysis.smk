"""Allele frequency analysis rules.

Two-stage flow:
1. ``compute_allele_freq_per_timepoint`` — runs once per (MAG, single timepoint)
   and writes a Parquet cache file with per-sample allele frequencies for that
   timepoint. Reused across every (timepoint_combination, group_combination)
   that includes that timepoint.
2. ``allele_analysis`` — runs once per (MAG, timepoint_combination, group_combination)
   and computes the per-combination diff/aggregate outputs from the two (or one,
   for single data) cache files.

The previously-emitted ``{mag}_allele_frequency_longitudinal.tsv.gz`` is no
longer produced; downstream rules read the Parquet cache files directly.
"""


def _allele_analysis_cache_input(wildcards):
    """Resolve the per-(gr_combo, timepoint) cache paths required for one combination."""
    if DATA_TYPE == "longitudinal":
        tps = wildcards.timepoints.split("_")
    else:
        tps = [wildcards.timepoints]
    return [
        get_allele_freq_cache_path(
            mag_wildcard=wildcards.mag,
            groups_wildcard=wildcards.groups,
            timepoint_wildcard=tp,
        )
        for tp in tps
    ]


rule compute_allele_freq_per_timepoint:
    input:
        # Single canonical QC file for this (gr_combo, timepoint).
        # Using ONE file (from the first tp_combo in config order that contains
        # this timepoint) means Snakemake resolves the dependency as a normal
        # qc + generate_metadata chain — no cross-combination checkpoint
        # dependencies that might not exist yet.
        qc_files=lambda wildcards: get_canonical_qc_file(
            wildcards.mag, wildcards.timepoint, wildcards.groups
        ),
    output:
        cache=os.path.join(
            OUTDIR,
            "allele_freq_cache",
            "{mag}_{groups}_{timepoint}_allele_frequency.parquet",
        ),
    params:
        data_type=DATA_TYPE,
        # For longitudinal data the script filters the QC file's rows to this
        # timepoint; for single data the QC file covers one timepoint already
        # so --timepoint is not passed.
        timepoint_arg=lambda wildcards: (
            f"--timepoint {wildcards.timepoint}"
            if DATA_TYPE == "longitudinal"
            else ""
        ),
    threads: get_threads("allele_freq_cache")
    retries: get_retries("allele_freq_cache")
    resources:
        mem_mb=get_mem_mb("allele_freq_cache"),
        time=get_time("allele_freq_cache"),
    shell:
        """
        alleleflux-cache-allele-freq \
            --magID {wildcards.mag} \
            --qc_files {input.qc_files} \
            {params.timepoint_arg} \
            --data_type {params.data_type} \
            --output_path {output.cache} \
            --cpus {threads}
        """


rule allele_analysis:
    input:
        cache_files=_allele_analysis_cache_input,
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
            ),
        ),
        allele_freq_no_zero_diff=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            (
                "{mag}_allele_frequency_no_constant.tsv.gz"
                if DATA_TYPE == "single"
                else "{mag}_allele_frequency_changes_no_zero-diff.tsv.gz"
            ),
        ) if not config["quality_control"].get("disable_zero_diff_filtering", False) else [],
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
        data_type=DATA_TYPE,
        groups_arg=lambda wildcards: "--groups " + " ".join(wildcards.groups.split("_")),
    threads: 1
    retries: get_retries("allele_analysis")
    resources:
        mem_mb=get_mem_mb("allele_analysis"),
        time=get_time("allele_analysis"),
    shell:
        """
        alleleflux-allele-freq \
            --magID {wildcards.mag} \
            --cache_files {input.cache_files} \
            --data_type {params.data_type} \
            --output_dir {params.outDir} \
            {params.groups_arg} \
            {params.disable_zero_diff_filtering}
        """
