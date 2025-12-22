"""
Sample metadata generation rules.

This module generates per-MAG metadata files that organize sample profile paths
by experimental group and timepoint. These metadata files are used as input for
QC analysis and downstream statistical tests.

When USE_EXISTING_PROFILES is True, profiles are read directly from PROFILES_DIR
(the existing profiles path specified in config). No profiling step runs.
The metadata files will contain the actual paths to the existing profile files.
"""

# Define input function to conditionally depend on profile rule outputs
def get_metadata_inputs(wildcards):
    """Return inputs for generate_metadata rule.
    
    When using existing profiles, we don't depend on the profile rule outputs.
    When generating new profiles, we depend on all sample directories.
    """
    inputs = {"metadata": config["input"]["metadata_path"]}
    
    if not USE_EXISTING_PROFILES:
        # Only depend on profile outputs when generating new profiles
        inputs["sampleDirs"] = expand(
            os.path.join(OUTDIR, "profiles", "{sample}"), sample=samples
        )
    
    return inputs


rule generate_metadata:
    input:
        unpack(get_metadata_inputs),
    output:
        outDir=directory(
            os.path.join(
                OUTDIR, "inputMetadata", "inputMetadata_{timepoints}-{groups}"
            )
        ),
    params:
        # Read profiles from PROFILES_DIR - either existing profiles or newly generated
        rootDir=PROFILES_DIR,
        data_type=DATA_TYPE,
        group_args=lambda wildcards: f"--groups {wildcards.groups.replace('_', ' ')}",
        timepoint_args=lambda wildcards: f"--timepoints {wildcards.timepoints.replace('_', ' ')}",
    resources:
        mem_mb=get_mem_mb("generate_metadata"),
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
