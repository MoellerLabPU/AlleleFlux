"""
Quality control rules.

This module performs QC on profiled samples, calculating coverage breadth and
depth metrics for each MAG. QC results determine sample eligibility for
downstream statistical analysis based on configurable thresholds.
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
        coverage_threshold=config["quality_control"]["coverage_threshold"],
        data_type=DATA_TYPE,
    threads: get_threads("qc")
    resources:
        time=get_time("qc"),
        mem_mb=get_mem_mb("qc"),
    shell:
        """
        alleleflux-qc \
            --rootDir {input.metadata_dir} \
            --fasta {input.fasta} \
            --breadth_threshold {params.breadth_threshold} \
            --coverage_threshold {params.coverage_threshold} \
            --cpus {threads} \
            --output_dir {output.outDir} \
            --data_type {params.data_type} \
            --mag_mapping_file {input.mag_mapping}
        """
