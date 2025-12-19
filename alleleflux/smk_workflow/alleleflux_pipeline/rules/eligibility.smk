"""
QC-based eligibility table checkpoint.

This checkpoint generates eligibility tables based on QC results, determining
which MAGs have sufficient sample coverage to proceed to downstream analysis.
The checkpoint mechanism allows Snakemake to dynamically update the DAG based
on which MAGs pass QC thresholds.
"""

checkpoint eligibility_table:
    input:
        qc_dir=os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}"),
    output:
        out_fPath=os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}.tsv"),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        data_type=DATA_TYPE,
    resources:
        time=get_time("eligibility_table"),
    shell:
        """
        alleleflux-eligibility \
            --qc_dir {input.qc_dir} \
            --min_sample_num {params.min_sample_num} \
            --output_file {output.out_fPath} \
            --data_type {params.data_type}
        """
