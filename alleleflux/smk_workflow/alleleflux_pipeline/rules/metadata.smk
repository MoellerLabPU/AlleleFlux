"""
Sample metadata generation rules.

This module generates per-MAG metadata files that organize sample profile paths
by experimental group and timepoint. These metadata files are used as input for
QC analysis and downstream statistical tests.

Note: When USE_EXISTING_PROFILES is True, profiles are read from PROFILES_DIR
(the existing profiles path), but the metadata is still written to OUTDIR.
"""

rule generate_metadata:
    input:
        metadata=config["input"]["metadata_path"],
        # This ensures that generate_metadata is run after all samples are profiled
        # Uses OUTDIR/profiles which will be symlinks when using existing profiles
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
        # Read profiles from the profiles directory in OUTDIR (which may be symlinks)
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
