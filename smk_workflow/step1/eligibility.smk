"""
Step 1 Eligibility Table Rules
Rules for generating eligibility tables for downstream analysis
"""

rule eligibility_table:
    input:
        qc_dir=os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}"),
    output:
        out_fPath=os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}.tsv"),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        data_type=DATA_TYPE,
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-eligibility \
            --qc_dir {input.qc_dir} \
            --min_sample_num {params.min_sample_num} \
            --output_file {output.out_fPath} \
            --data_type {params.data_type}
        """
