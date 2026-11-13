"""
Quality control rules.

This module performs QC on profiled samples, calculating coverage breadth and
depth metrics for each MAG. QC results determine sample eligibility for
downstream statistical analysis based on configurable thresholds.

**Scope: one job per single timepoint, every group.**
QC metrics (breadth, depth) are intrinsic to each sample and independent of
which comparison or group pair will later use them.  The QC rule runs once per
{timepoint} wildcard (a single timepoint), on the metadata folder of that
timepoint, and writes PER-SAMPLE results only: one ``{mag}_QC.tsv`` per genome
of the run, header-only when the genome has no sample at this timepoint.

Which subjects are readable at both timepoints of a comparison, and how many
replicates each tested group then has, is computed downstream in the
eligibility_table checkpoint from the QC folders of the comparison's two
timepoints.
"""

rule qc:
    """Run QC once per single timepoint, covering ALL groups.

    Input is that timepoint's metadata directory (no {groups} wildcard).
    Output is that timepoint's QC directory, reused by every comparison naming
    the timepoint (eligibility_table checkpoints, allele-frequency caches,
    pairwise ANI).
    """
    input:
        # Depend on the metadata sentinel rather than the directory itself.
        metadata_done=metadata_sentinel("{timepoint}"),
        fasta=config["input"]["fasta_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        # Sentinel marker (see common.smk) — replaces directory() output.
        sentinel=qc_sentinel("{timepoint}"),
    params:
        breadth_threshold=config["quality_control"]["breadth_threshold"],
        coverage_threshold=config["quality_control"]["coverage_threshold"],
        data_type=DATA_TYPE,
        # Concrete directory paths the script reads from / writes into.
        metadata_dir=os.path.join(
            OUTDIR, "inputMetadata", "inputMetadata_{timepoint}"
        ),
        outDir=os.path.join(OUTDIR, "QC", "QC_{timepoint}"),
    threads: get_threads("qc")
    retries: get_retries("qc")
    resources:
        time=get_time("qc"),
        runtime=get_runtime("qc"),
        mem_mb=get_mem_mb("qc"),
    shell:
        """
        alleleflux-qc \
            --rootDir {params.metadata_dir} \
            --fasta {input.fasta} \
            --breadth_threshold {params.breadth_threshold} \
            --coverage_threshold {params.coverage_threshold} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --data_type {params.data_type} \
            --mag_mapping_file {input.mag_mapping}
        touch {output.sentinel}
        """
