"""
Step 1 Metadata Generation Rules
Rules for generating metadata from profiled samples
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
        time=config["resources"]["time"]["general"],
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
