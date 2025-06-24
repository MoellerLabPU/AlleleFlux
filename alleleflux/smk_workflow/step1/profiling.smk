"""
Step 1 Profiling Rules
Rules for sample profiling based on BAM files
"""

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
