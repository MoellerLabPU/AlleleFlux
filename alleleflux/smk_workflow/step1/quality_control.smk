"""
Step 1 Quality Control Rules
Rules for quality control of profiled samples
"""

rule qc:
    input:
        metadata_dir=os.path.join(
            OUTDIR, "inputMetadata", "inputMetadata_{timepoints}-{groups}"
        ),
        fasta=config["input"]["fasta_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        outDir=directory(os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}")),
    params:
        breadth_threshold=config["quality_control"]["breadth_threshold"],
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["threads_per_job"]
    resources:
        time=config["resources"]["time"]["general"],
        mem_mb=config["resources"]["memory"]["quality_control"],
    shell:
        """
        alleleflux-qc \
            --rootDir {input.metadata_dir} \
            --fasta {input.fasta} \
            --breadth_threshold {params.breadth_threshold} \
            --cpus {threads} \
            --output_dir {output.outDir} \
            --data_type {params.data_type} \
            --mag_mapping_file {input.mag_mapping}
        """
