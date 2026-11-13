"""
Sample metadata generation rules.

This module generates per-MAG metadata files that organize sample profile paths
by experimental group and timepoint. These metadata files are used as input for
QC analysis and downstream statistical tests.

**Scope: one job per single timepoint, every group.**
The metadata and QC rules run once per single timepoint (``pre``), not per
comparison (``pre_post``) and not per group pair, passing ALL unique group
values from the config.  A timepoint shared by several comparisons is read
once: comparisons pre_post and pre_end need 3 metadata + 3 QC jobs (pre, post,
end), where one job per comparison would read the ``pre`` samples twice.

Every genome with a profile anywhere in the run gets a ``{mag}_metadata.tsv``
in every timepoint's folder, header-only when it has no sample there.

Group-pair-specific eligibility (which MAGs have enough samples in group A vs B)
is computed per (timepoints, groups) in the eligibility_table checkpoint, which
reads the QC folders of the comparison's two timepoints.

When USE_EXISTING_PROFILES is True, profiles are read directly from PROFILES_DIR
(the existing profiles path specified in config). No profiling step runs.
The metadata files will contain the actual paths to the existing profile files.
"""

# ALL unique group values across all group combinations — passed to metadata
# so a single job covers every group in this timepoint, not just one pair.
_all_group_values = sorted(
    {gr["treatment"] for gr in config["analysis"]["groups_combinations"]}
    | {gr["control"] for gr in config["analysis"]["groups_combinations"]}
)


# Define input function to conditionally depend on profile rule outputs
def get_metadata_inputs(wildcards):
    """Return inputs for generate_metadata rule.
    
    When using existing profiles, we don't depend on the profile rule outputs.
    When generating new profiles, we depend on all sample directories.
    """
    inputs = {"metadata": config["input"]["metadata_path"]}
    
    if not USE_EXISTING_PROFILES:
        # Only depend on profile outputs when generating new profiles.
        # Depend on the per-sample profile sentinel — see common.smk for the
        # rationale.  Using the sentinel (rather than directory()) means an
        # interrupted profile job cannot silently invalidate downstream rules.
        inputs["sampleDirs"] = [profile_sentinel(s) for s in samples]

    return inputs


rule generate_metadata:
    """Generate per-MAG metadata files for ONE timepoint.

    Runs once per {timepoint} wildcard (a single timepoint, not a comparison
    and not a group pair).  All unique group values from the config are passed
    via --groups so downstream QC sees every group in one directory.
    """
    input:
        unpack(get_metadata_inputs),
    output:
        # Sentinel marker (see common.smk) — replaces a directory() output to
        # avoid pre-deletion of the inputMetadata dir on every re-run.
        sentinel=metadata_sentinel("{timepoint}"),
    params:
        # Read profiles from PROFILES_DIR - either existing profiles or newly generated
        rootDir=PROFILES_DIR,
        data_type=DATA_TYPE,
        # The actual output directory the script writes per-MAG metadata into.
        outDir=os.path.join(
            OUTDIR, "inputMetadata", "inputMetadata_{timepoint}"
        ),
        # Pass every unique group so the QC directory covers all groups
        group_args=" ".join(_all_group_values),
        # Exactly one timepoint per job: the samples of this timepoint only.
        timepoint_args=lambda wildcards: f"--timepoints {wildcards.timepoint}",
    retries: get_retries("generate_metadata")
    resources:
        mem_mb=get_mem_mb("generate_metadata"),
        time=get_time("generate_metadata"),
        runtime=get_runtime("generate_metadata"),
    shell:
        """
        alleleflux-metadata \
            --rootDir {params.rootDir} \
            --metadata {input.metadata} \
            --outDir {params.outDir} \
            --groups {params.group_args} \
            --data_type {params.data_type} \
            {params.timepoint_args}
        touch {output.sentinel}
        """
