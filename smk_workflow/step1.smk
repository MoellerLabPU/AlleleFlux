"""
Workflow for Step 1: Sample Profiling and Eligibility Table Generation
This workflow is designed to profile samples based on BAM files and generate eligibility table
to be used by step 2.
"""


# Include modular workflow components
include: "rules/common.smk"


samples = [
    os.path.basename(bam).split(".")[0]
    for bam in glob(os.path.join(config["bamDir"], "*.sorted.bam"))
]


localrules:
    all,


rule all:
    input:
        expand(os.path.join(OUTDIR, "profiles", "{sample}.sorted"), sample=samples),
        expand(
            os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}.tsv"),
            timepoints=timepoints_labels,
            groups=groups_labels,
        ),


#####################################
# Step 1: Profile Samples
#####################################


rule profile:
    input:
        bam=os.path.join(config["bamDir"], "{sample}.sorted.bam"),
        fasta=config["fasta"],
        prodigal=config["prodigal"],
    output:
        sampleDirs=directory(os.path.join(OUTDIR, "profiles", "{sample}.sorted")),
    threads: config["cpus"]["profile"]
    resources:
        mem_mb=config["memory"]["profile"],
        time=config["time"]["profile"],
    params:
        outDir=os.path.join(OUTDIR, "profiles"),
        scriptPath=config["scripts"]["profile"],
    benchmark:
        os.path.join(
            OUTDIR,
            "benchmarks",
            "profile_{sample}.tsv",
        )
    shell:
        """
        python {params.scriptPath} \
        --bam_path {input.bam} --fasta_path {input.fasta} --prodigal_fasta {input.prodigal} \
        --cpus {threads} --output_dir {params.outDir}
        """


rule generate_metadata:
    input:
        metadata=config["metadata_file"],
        # This enures that generate_metadata is run after all samples are profiled
        sampleDirs=expand(
            os.path.join(OUTDIR, "profiles", "{sample}.sorted"), sample=samples
        ),
    output:
        outDir=directory(
            os.path.join(
                OUTDIR, "inputMetadata", "inputMetadata_{timepoints}-{groups}"
            )
        ),
    params:
        rootDir=os.path.join(OUTDIR, "profiles"),
        scriptPath=config["scripts"]["generate_mag_metadata"],
        data_type=DATA_TYPE,
        group_args=lambda wildcards: f"--groups {wildcards.groups.replace('_', ' ')}",
        timepoint_args=lambda wildcards: f"--timepoints {wildcards.timepoints.replace('_', ' ')}",
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
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
        fasta=config["fasta"],
    output:
        outDir=directory(os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}")),
    params:
        breadth_threshold=config["breadth_threshold"],
        scriptPath=config["scripts"]["quality_control"],
        data_type=DATA_TYPE,
    threads: config["cpus"]["quality_control"]
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
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
        min_sample_num=config["min_sample_num"],
        script=config["scripts"]["eligibility_table"],
        data_type=DATA_TYPE,
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.script} \
            --qc_dir {input.qc_dir} \
            --min_sample_num {params.min_sample_num} \
            --output_file {output.out_fPath} \
            --data_type {params.data_type}
        """
