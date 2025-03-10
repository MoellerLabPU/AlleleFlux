import os
import pandas as pd
from glob import glob
from collections import defaultdict
from multiprocessing import Pool
import logging


# Load the configuration file
configfile: os.path.join(workflow.basedir, "config.yml")


OUTDIR = config["root_out"]
# Automatically detect sample IDs based on filenames in the directory
samples = [
    os.path.basename(bam).split(".")[0]
    for bam in glob(os.path.join(config["bamDir"], "*.sorted.bam"))
]


timepoints_labels = [
    f"{tp[0]}_{tp[1]}" for tp in config["timepoints_combinations"]
]  # ["T1_T2", "T2_T3"]
groups_labels = [
    f"{gr[0]}_{gr[1]}" for gr in config["groups_combinations"]
]  # ["G1_G2", "G2_G4"]


wildcard_constraints:
    groups="[a-zA-Z0-9_-]+",  # Allow alphanumeric, underscores, and hyphens
    timepoints="[a-zA-Z0-9_-]+",


localrules:
    all,


rule all:
    input:
        expand(os.path.join(OUTDIR, "profiles", "{sample}.sorted"), sample=samples),
        expand(
            os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}"),
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
    resources:
        time=config["time"]["general"],
    # log:
    #     os.path.join(config["logDir"], "generate_metadata", "{timepoints}-{groups}.log"),
    shell:
        """
        # Extract the two timepoints and groups from the wildcard strings
        TP1=$(echo {wildcards.timepoints} | cut -d'_' -f1)
        TP2=$(echo {wildcards.timepoints} | cut -d'_' -f2)
        G1=$(echo {wildcards.groups} | cut -d'_' -f1)
        G2=$(echo {wildcards.groups} | cut -d'_' -f2)

        python {params.scriptPath} \
            --rootDir {params.rootDir} \
            --metadata {input.metadata} \
            --outDir {output.outDir} \
            --timepoints $TP1 $TP2 \
            --groups $G1 $G2 
         
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
            --output_dir {output.outDir}
        """


rule eligibility_table:
    input:
        qc_dir=os.path.join(OUTDIR, "QC", "QC_{timepoints}-{groups}"),
    output:
        out_fPath=os.path.join(OUTDIR, "eligibility_table_{timepoints}-{groups}"),
    params:
        min_sample_num=config["min_sample_num"],
        script=config["scripts"]["eligibility_table"],
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.script} \
            --qc_dir {input.qc_dir} \
            --min_sample_num {params.min_sample_num} \
            --output_file {output.out_fPath}
        """
