"""
QC-based eligibility table checkpoint.

This checkpoint generates eligibility tables based on QC results, determining
which MAGs have sufficient sample coverage to proceed to downstream analysis.
The checkpoint mechanism allows Snakemake to dynamically update the DAG based
on which MAGs pass QC thresholds.

**Per-timepoint QC, per-comparison eligibility:**
QC runs once per single timepoint and holds per-sample results for every
group.  The eligibility_table checkpoint runs once per (timepoints, groups)
pair: it reads the QC folders of the comparison's timepoints (two for
longitudinal data, one for single), works out which subjects are readable at
both timepoints, and counts replicates over the two tested groups given by
--groups.
"""

checkpoint eligibility_table:
    input:
        # Depend on the QC sentinel of EACH timepoint of the comparison (two
        # for longitudinal data, one for single).  The QC directories are
        # consumed by the script, but the sentinels are what prove QC finished
        # — see common.smk for the rationale.
        qc_done=lambda wildcards: [
            qc_sentinel(tp) for tp in comparison_timepoints(wildcards.timepoints)
        ],
        # Permuted (null) run: depend on the per-pair sheet so it is built before
        # this checkpoint (generate mode only; [] otherwise — no behaviour change).
        permuted_md=lambda wildcards: permuted_metadata_input(wildcards.groups),
    output:
        out_fPath=os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}.tsv"),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        data_type=DATA_TYPE,
        # The two tested groups: replicates are counted over these two only.
        groups_args=lambda wildcards: wildcards.groups.replace("_", " "),
        # One QC folder per timepoint of the comparison, earlier first.
        # qc_dir_of resolves under REUSE_DIR (== OUTDIR normally); for a
        # reuse_from (permuted) run that is the real run's QC, so this
        # checkpoint reads the existing QC.
        qc_dirs=lambda wildcards: " ".join(
            qc_dir_of(tp) for tp in comparison_timepoints(wildcards.timepoints)
        ),
        # Permuted (null) run: relabel the groups of the reused QC rows from the
        # permuted sheet before the replicates are counted.
        permuted_metadata=lambda wildcards: permuted_metadata_flag(wildcards.groups),
    retries: get_retries("eligibility_table")
    resources:
        time=get_time("eligibility_table"),
        runtime=get_runtime("eligibility_table"),
    shell:
        """
        alleleflux-eligibility \
            --qc_dirs {params.qc_dirs} \
            --min_sample_num {params.min_sample_num} \
            --output_file {output.out_fPath} \
            --data_type {params.data_type} \
            --groups {params.groups_args} \
            {params.permuted_metadata}
        """
