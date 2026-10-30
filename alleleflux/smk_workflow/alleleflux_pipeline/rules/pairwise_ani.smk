"""Pairwise conANI / popANI between the QC-passing sample pairs of each MAG.

Scope: ONE job per {mag}, independent of group and timepoint combinations.  The
sample universe is the union of everything that passed QC at any timepoint
combination, so a single job covers within-mouse (pre vs end) and between-mouse
comparisons at once; which pairs are actually compared is the CLI's --pairs
setting (config analysis.pairwise_ani.pairs).

The MAG universe is get_tested_mags(): every MAG eligible for at least one
enabled test in at least one comparison, read from the eligibility checkpoints
(so the targets are generated in get_final_pipeline_outputs, after the loop that
forces those checkpoints).  Untested MAGs never get an ANI job.  The rule's own
inputs are the per-timepoint QC sentinels plus the MAG's allele-frequency cache
files (the single per-site source after QC), so the order is QC -> cache -> ANI.
"""


rule pairwise_ani:
    input:
        # Depend on the QC sentinels rather than the QC TSVs so an interrupted QC
        # job cannot leave a half-written directory that looks complete.
        qc_sentinels=[qc_sentinel(tp) for tp in timepoints_labels],
        # The MAG's cache files (one per timepoint it is QC-eligible at): the
        # single per-site source after QC, already site-depth-filtered.
        caches=lambda wildcards: allele_freq_cache_files_for_mag(wildcards.mag),
        fasta=config["input"]["fasta_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        # All three files the CLI writes, so Snakemake tracks (and cleans up on
        # failure) the complete set, not just the primary table.  The primary
        # path comes from the shared helper (called bare, it returns the "{mag}"
        # wildcard pattern) so the rule, the target generator, and any future
        # consumer rule can never spell it differently.
        pair_table=get_pairwise_ani_output_path(),
        sample_table=os.path.join(
            OUTDIR, "pairwise_ani", "{mag}_pairwise_ani_samples.tsv"
        ),
    retries: get_retries("pairwise_ani")
    # Single-threaded end to end: the cache is streamed serially (fast, and a
    # pool would double peak memory) and the comparison loop is single-threaded.
    threads: 1
    resources:
        mem_mb=get_mem_mb("pairwise_ani"),
        time=get_time("pairwise_ani"),  # HH:MM:SS -- cluster-generic profile
        runtime=get_runtime("pairwise_ani"),  # minutes  -- native slurm plugin
    params:
        # Wildcard-dependent params need a lambda (bare {wildcards.X} does not
        # interpolate inside params values -- house rule).
        qc_files=lambda wildcards: " ".join(get_all_qc_files_for_mag(wildcards.mag)),
        cache_dir=os.path.join(REUSE_DIR, "allele_freq_cache"),
        output_dir=os.path.join(OUTDIR, "pairwise_ani"),
        # The presence rule, shared with baseline_presence via one dict in
        # common.smk so the two can never disagree about what "present" means.
        min_cov=PRESENCE_RULE["min_cov"],
        min_freq=PRESENCE_RULE["min_freq"],
        fdr=PRESENCE_RULE["fdr"],
        min_base_quality=PRESENCE_RULE["min_base_quality"],
        # The two mode knobs come from PAIRWISE_ANI (common.smk), resolved once
        # with data-type-aware defaults: single data gets pairs=all and no SNP
        # locations, and the same-mouse modes are refused there at parse time.
        # PAIRWISE_ANI is None when the feature is off; the ``or {}`` keeps this
        # params block parseable then (no job of this rule exists in that case).
        pairs=(PAIRWISE_ANI or {}).get("pairs", ""),
        # --transitions is only meaningful for pairs == "transitions"; built from the
        # configured timepoint combinations (EARLIER:LATER each).  The guard in
        # common.smk refuses "transitions" on single data, so every entry here
        # has two timepoints when this branch is taken.
        transitions_arg=(
            "--transitions " + " ".join(
                f"{tc['timepoint'][0]}:{tc['timepoint'][1]}"
                for tc in config["analysis"]["timepoints_combinations"]
            )
            if (PAIRWISE_ANI or {}).get("pairs") == "transitions"
            else ""
        ),
        store_snp_locations=(PAIRWISE_ANI or {}).get("store_snp_locations", ""),
    shell:
        """
        alleleflux-pairwise-ani \
            --mag {wildcards.mag} \
            --cache_dir {params.cache_dir} \
            --qc_files {params.qc_files} \
            --fasta {input.fasta} \
            --mag_mapping {input.mag_mapping} \
            --output_dir {params.output_dir} \
            --min_cov {params.min_cov} \
            --min_freq {params.min_freq} \
            --fdr {params.fdr} \
            --min_base_quality {params.min_base_quality} \
            --pairs {params.pairs} \
            {params.transitions_arg} \
            --store_snp_locations {params.store_snp_locations}
        """
