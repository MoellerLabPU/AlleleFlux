"""Was each significant allele already present at the baseline timepoint?

``baseline_presence`` runs ONCE per {timepoints}-{groups} comparison, on the
p_value_summary file of the configured test family (analysis.baseline_presence.summary)
and the per-MAG significance-test files behind it, and reads every sample's
profile once (both groups, both timepoints).  It therefore sits AFTER the
statistics stage: its targets are generated inside get_final_pipeline_outputs
(checkpoint-aware), only where the family's summary target itself exists.

Writes three files: the long table, the summary, and the rising table (only
the allele that rose, per group, judged on that group's later samples).

When use_strain_turnover is on, the turnover tables are joined in for the
strain_background column and become an input dependency.
"""


def _baseline_presence_turnover_inputs(wildcards):
    """Turnover tables as inputs only when the strain column is requested."""
    if config["analysis"].get("use_strain_turnover", False):
        return generate_strain_turnover_targets()
    return []


rule baseline_presence:
    input:
        summary=os.path.join(
            OUTDIR,
            "p_value_summary",
            "{timepoints}-{groups}",
            "p_value_summary_" + BASELINE_PRESENCE_FAMILY + "_{timepoints}-{groups}.tsv",
        ),
        turnover=_baseline_presence_turnover_inputs,
        fasta=config["input"]["fasta_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
        metadata=config["input"]["metadata_path"],
    output:
        long=get_baseline_presence_stem() + ".tsv.gz",
        summary=get_baseline_presence_stem() + "_summary.tsv",
        rising=get_baseline_presence_stem() + "_rising.tsv",
    retries: get_retries("baseline_presence")
    threads: get_threads("baseline_presence")
    resources:
        mem_mb=get_mem_mb("baseline_presence"),
        time=get_time("baseline_presence"),
        runtime=get_runtime("baseline_presence"),
    params:
        run_dir=OUTDIR,
        profiles_dir=PROFILES_DIR,
        output_dir=os.path.join(OUTDIR, "baseline_presence"),
        summary_family=BASELINE_PRESENCE_FAMILY,
        # the shared constant (common.smk), not a direct config read: params are
        # evaluated for every config, including ones without a baseline_presence block
        test_type=BASELINE_PRESENCE_TEST_TYPE,
        threshold_column=config["analysis"].get("baseline_presence", {}).get("threshold_column", "q_value"),
        threshold=config["analysis"].get("baseline_presence", {}).get("threshold", 0.05),
        # which samples a site is judged on (comparison | own_group); the CLI's own
        # default when the config does not say.  own_group needs a within-group
        # summary family (single_sample / lmm_across_time); the CLI refuses the rest.
        baseline_scope=config["analysis"].get("baseline_presence", {}).get("baseline_scope", "comparison"),
        # Per-site depth filter flags, shared with allele_freq_cache and
        # pairwise_ani so every profile reader filters identically.
        site_filter_args=SITE_DEPTH_FILTER_ARGS,
        # The presence rule, the SAME dict pairwise_ani reads (common.smk), so
        # "present" means one thing across the ANI table and this one.
        min_cov=PRESENCE_RULE["min_cov"],
        min_freq=PRESENCE_RULE["min_freq"],
        fdr=PRESENCE_RULE["fdr"],
        min_base_quality=PRESENCE_RULE["min_base_quality"],
        turnover_arg=(
            f"--turnover_dir {os.path.join(OUTDIR, 'strain_turnover')}"
            if config["analysis"].get("use_strain_turnover", False)
            else ""
        ),
    shell:
        """
        alleleflux-baseline-presence \
            --run_dir {params.run_dir} \
            --comparison {wildcards.timepoints}-{wildcards.groups} \
            --summary {params.summary_family} \
            --test_type {params.test_type} \
            --threshold_column {params.threshold_column} \
            --threshold {params.threshold} \
            --baseline_scope {params.baseline_scope} \
            --profiles_dir {params.profiles_dir} \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --mag_mapping {input.mag_mapping} \
            --output_dir {params.output_dir} \
            --min_cov {params.min_cov} \
            --min_freq {params.min_freq} \
            --fdr {params.fdr} \
            --min_base_quality {params.min_base_quality} \
            {params.turnover_arg} \
            {params.site_filter_args} \
            --cpus {threads}
        """


def _rising_origin_inputs(wildcards):
    """Every comparison's rising table, plus the classification when strain turnover runs.

    Reads the same list generate_rising_origin_targets gates on, so the rule is
    only ever scheduled with at least one rising table to read.
    """
    return generate_rising_table_targets() + generate_replacement_classification_targets()


# Seconds of pandas over small tables: not worth a cluster job.
localrules:
    rising_origin,


rule rising_origin:
    """Standing vs de novo for the alleles that rose, over every comparison: one
    tidy summary.  With strain turnover on, de novo in a replaced MAG becomes
    new_strain (the run's classification table, its configured voter floor)."""
    input:
        _rising_origin_inputs,
    output:
        alleles=os.path.join(OUTDIR, "baseline_presence", "rising_origin", "rising_origin_alleles.tsv"),
        summary=os.path.join(OUTDIR, "baseline_presence", "rising_origin", "rising_origin_summary.tsv"),
    params:
        rising_dir=os.path.join(OUTDIR, "baseline_presence"),
        output_dir=os.path.join(OUTDIR, "baseline_presence", "rising_origin"),
        classification_arg=(
            f"--classification {get_replacement_classification_path()}"
            if config["analysis"].get("use_strain_turnover", False)
            else ""
        ),
    shell:
        """
        alleleflux-rising-origin \
            --rising_dir {params.rising_dir} \
            --output_dir {params.output_dir} \
            {params.classification_arg}
        """
