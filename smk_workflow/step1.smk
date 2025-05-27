"""
Workflow for Step 1: Sample Profiling and Eligibility Table Generation
This workflow is designed to profile samples based on BAM files and generate eligibility table
to be used by step 2.
"""

include: "shared/common.smk"
# Get sample information from metadata file
samples, sample_to_bam_map = get_sample_info()

# Include modular workflow components
include: "step1/profiling.smk"
include: "step1/metadata.smk"
include: "step1/quality_control.smk"
include: "step1/eligibility.smk"


localrules:
    all,

rule all:
    input:
        expand(os.path.join(OUTDIR, "profiles", "{sample}"), sample=samples),
        expand(
            os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}.tsv"),
            timepoints=timepoints_labels,
            groups=groups_labels,
        ),
