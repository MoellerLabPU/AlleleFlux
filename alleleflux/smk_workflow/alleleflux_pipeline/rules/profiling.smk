"""
Sample profiling rules.

This module contains rules for profiling BAM files to extract allele frequency
information at each genomic position. The profiling step is the foundation of
the AlleleFlux pipeline.

Reference files (FASTA, Prodigal, MAG mapping) are wrapped with ancient() to
prevent unnecessary re-runs when these stable files have updated timestamps.

If USE_EXISTING_PROFILES is True (profiles_path specified in config), the
profiling step is skipped and existing profiles are used instead. This enables
reusing profiles across multiple analysis runs with different parameters.

Expected directory structure for profiles_path:
    {profiles_path}/
        {sample1}/
            {sample1}_{MAG1}_profiled.tsv.gz
            {sample1}_{MAG2}_profiled.tsv.gz
            ...
        {sample2}/
            {sample2}_{MAG1}_profiled.tsv.gz
            ...
"""

if USE_EXISTING_PROFILES:
    # When using existing profiles, create a validation rule that checks
    # the profile directory exists and creates a symlink in the output directory
    rule use_existing_profiles:
        """
        Validate and link existing profile directories.
        
        This rule is used when profiles_path is specified in the config.
        It validates that the expected sample subdirectory exists in the
        pre-existing profiles location (profiles_path/{sample}/) and creates
        a symbolic link to it in the output directory.
        
        Expected structure in profiles_path:
            profiles_path/{sample}/{sample}_{MAG}_profiled.tsv.gz
        """
        input:
            existing_profile=lambda wildcards: os.path.join(EXISTING_PROFILES_PATH, wildcards.sample),
        output:
            sampleDirs=directory(os.path.join(OUTDIR, "profiles", "{sample}")),
        run:
            import os
            import shutil
            
            src_dir = input.existing_profile
            dst_dir = output.sampleDirs
            
            # Validate the source directory exists
            if not os.path.isdir(src_dir):
                raise FileNotFoundError(
                    f"Profile directory not found: {src_dir}. "
                    f"Ensure profiles_path in config is correct and contains "
                    f"profiles for sample '{wildcards.sample}'."
                )
            
            # Check that the directory contains profile files
            profile_files = list(Path(src_dir).glob("*_profiled.tsv.gz"))
            if not profile_files:
                raise FileNotFoundError(
                    f"No profile files (*_profiled.tsv.gz) found in {src_dir}. "
                    f"The directory may be empty or contain invalid profiles."
                )
            
            # Create parent directory if needed
            os.makedirs(os.path.dirname(dst_dir), exist_ok=True)
            
            # Create symbolic link to the existing profile directory
            # This avoids copying data while making it accessible at the expected path
            if os.path.islink(dst_dir):
                os.unlink(dst_dir)
            elif os.path.exists(dst_dir):
                shutil.rmtree(dst_dir)
            
            os.symlink(os.path.abspath(src_dir), dst_dir)

else:
    # Standard profiling rule - generates new profiles from BAM files
    rule profile:
        input:
            bam=lambda wildcards: sample_to_bam_map[wildcards.sample],
            fasta=config["input"]["fasta_path"],
            prodigal=config["input"]["prodigal_path"],
            mag_mapping=config["input"]["mag_mapping_path"],
        output:
            sampleDirs=directory(os.path.join(OUTDIR, "profiles", "{sample}")),
        threads: get_threads("profile")
        resources:
            mem_mb=get_mem_mb("profile"),
            time=get_time("profile"),
        params:
            outDir=os.path.join(OUTDIR, "profiles"),
            no_ignore_orphans="--no-ignore-orphans" if not config.get("profiling", {}).get("ignore_orphans", True) else "",
            no_ignore_overlaps="--no-ignore-overlaps" if not config.get("profiling", {}).get("ignore_overlaps", True) else "",
            min_base_quality=config.get("profiling", {}).get("min_base_quality", 30),
            min_mapping_quality=config.get("profiling", {}).get("min_mapping_quality", 2),
            log_level=config.get("log_level", "INFO"),
        shell:
            """
            alleleflux-profile \
                --bam_path {input.bam} \
                --fasta_path {input.fasta} \
                --prodigal_fasta {input.prodigal} \
                --mag_mapping_file {input.mag_mapping} \
                --cpus {threads} \
                --output_dir {params.outDir} \
                --sampleID {wildcards.sample} \
                {params.no_ignore_orphans} \
                {params.no_ignore_overlaps} \
                --min-base-quality {params.min_base_quality} \
                --min-mapping-quality {params.min_mapping_quality} \
                --log-level {params.log_level}
            """
