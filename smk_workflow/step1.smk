"""
Workflow for Step 1: Sample Profiling and Eligibility Table Generation
This workflow is designed to profile samples based on BAM files and generate eligibility table
to be used by step 2.
"""

# Include modular workflow components
include: "rules/common.smk"

# Get sample information from metadata file
samples, sample_to_bam_map = get_sample_info()

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

rule profile:
    input:
        bam=lambda wildcards: sample_to_bam_map[wildcards.sample],
        fasta=config["input"]["fasta_path"],
        prodigal=config["input"]["prodigal_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        sampleDirs=directory(os.path.join(OUTDIR, "profiles", "{sample}")),
    threads: config["resources"]["cpus"]["profile"]
    resources:
        mem_mb=config["resources"]["memory"]["profile"],
        time=config["resources"]["time"]["profile"],
    params:
        outDir=os.path.join(OUTDIR, "profiles"),
    shell:
        """
        alleleflux-profile \
            --bam_path {input.bam} \
            --fasta_path {input.fasta} \
            --prodigal_fasta {input.prodigal} \
            --mag_mapping_file {input.mag_mapping} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            --sampleID {wildcards.sample}
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

rule qc:
    input:
        metadata_dir=os.path.join(
            OUTDIR, "inputMetadata", "inputMetadata_{timepoints}-{groups}"
        ),
        fasta=config["input"]["fasta_path"],
    output:
        outDir=directory(os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}")),
    params:
        breadth_threshold=config["quality_control"]["breadth_threshold"],
        data_type=DATA_TYPE,
    threads: config["resources"]["cpus"]["quality_control"]
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-qc \
            --rootDir {input.metadata_dir} \
            --fasta {input.fasta} \
            --breadth_threshold {params.breadth_threshold} \
            --cpus {threads} \
            --output_dir {output.outDir} \
            --data_type {params.data_type}
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
