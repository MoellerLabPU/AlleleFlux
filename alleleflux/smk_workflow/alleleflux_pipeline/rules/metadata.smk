"""
Sample metadata generation rules.

This module generates per-MAG metadata files that organize sample profile paths
by experimental group and timepoint. These metadata files are used as input for
QC analysis and downstream statistical tests.
"""

rule generate_metadata:
    input:
        metadata=config["input"]["metadata_path"],
        # This ensures that generate_metadata is run after all samples are profiled
        sampleDirs=expand(
            os.path.join(OUTDIR, "profiles", "{sample}"), sample=samples
        ),
    output:
        outDir=directory(
            os.path.join(
                OUTDIR, "inputMetadata", "inputMetadata_{timepoints}-{groups}"
            )
        ),
    params:
        rootDir=os.path.join(OUTDIR, "profiles"),
        data_type=DATA_TYPE,
        group_args=lambda wildcards: f"--groups {wildcards.groups.replace('_', ' ')}",
        timepoint_args=lambda wildcards: f"--timepoints {wildcards.timepoints.replace('_', ' ')}",
    resources:
        time=get_time("generate_metadata"),
    shell:
        """
        alleleflux-metadata \
            --rootDir {params.rootDir} \
            --metadata {input.metadata} \
            --outDir {output.outDir} \
            {params.group_args} \
            --data_type {params.data_type} \
            {params.timepoint_args}
        """
