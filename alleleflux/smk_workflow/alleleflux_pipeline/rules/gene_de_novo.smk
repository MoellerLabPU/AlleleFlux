"""Gene-level de novo accumulation: which genes collect new alleles, per cage.

``gene_de_novo`` runs ONCE per {timepoints}-{groups} comparison.  It reads the
comparison's mean-change tables (one per QC-eligible MAG: the candidate
rises), the allele-frequency cache (readable bases per cage from the paired
mice; the verdicts of every loaded sample at the candidate sites), the gene
file and, under the carrier strain gate, the per-MAG strain-turnover tables.
It writes three files per MAG (events, gene x cage, cage totals) next to one
tests table per comparison (one row per MAG x gene: Test A per group, the
gene against its own genome; Test B between the groups; BH q-values).

Its targets are generated with the other downstream analyses (gated by
use_gene_de_novo, off under allele_analysis_only, longitudinal only, at
least one QC-eligible MAG).  The knobs are read once in common.smk
(``gene_de_novo_params``); the cage floor is quality_control.min_sample_num
and BH follows statistics.fdr_group_by_mag_id, both shared with the rest of
the run.
"""


def _gene_de_novo_turnover_inputs(wildcards):
    """The per-MAG turnover tables as inputs under the carrier gate only."""
    if (GENE_DE_NOVO_PARAMS or {}).get("strain_gate") == "carrier":
        return generate_strain_turnover_targets()
    return []


rule gene_de_novo:
    input:
        # The candidate rises: every QC-eligible MAG's mean-change table of
        # this comparison (the rule reads the folder; the targets make the
        # dependency explicit so a stale table re-runs the allele step).
        mean_changes=lambda wildcards: generate_allele_analysis_targets(
            wildcards.timepoints, wildcards.groups
        ),
        # The comparison's cache files (both timepoints, every eligible MAG).
        caches=lambda wildcards: allele_freq_cache_files_for_comparison(
            wildcards.timepoints, wildcards.groups
        ),
        turnover=_gene_de_novo_turnover_inputs,
        metadata=config["input"]["metadata_path"],
        prodigal_fasta=config["input"]["prodigal_path"],
        fasta=config["input"]["fasta_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        tests=get_gene_de_novo_tests_path(),
    retries: get_retries("gene_de_novo")
    threads: get_threads("gene_de_novo")
    resources:
        mem_mb=get_mem_mb("gene_de_novo"),
        time=get_time("gene_de_novo"),
        runtime=get_runtime("gene_de_novo"),
    params:
        cache_dir=os.path.join(REUSE_DIR, "allele_freq_cache"),
        allele_analysis_dir=lambda wildcards: os.path.join(
            OUTDIR,
            "allele_analysis",
            f"allele_analysis_{wildcards.timepoints}-{wildcards.groups}",
        ),
        output_dir=lambda wildcards: os.path.join(
            OUTDIR, "gene_de_novo", f"{wildcards.timepoints}-{wildcards.groups}"
        ),
        # the shared de novo rules + this analysis's own knobs, one string (common.smk)
        knob_args=GENE_DE_NOVO_ARGS,
        # the cages-per-group floor is the run's replicate floor (decision 8)
        min_cages=config["quality_control"]["min_sample_num"],
        fdr_arg=(
            "--fdr_group_by_mag_id"
            if config["statistics"].get("fdr_group_by_mag_id", False)
            else ""
        ),
        # the presence rule, the SAME dict pairwise_ani and baseline_presence read
        min_cov=PRESENCE_RULE["min_cov"],
        min_freq=PRESENCE_RULE["min_freq"],
        fdr=PRESENCE_RULE["fdr"],
        min_base_quality=PRESENCE_RULE["min_base_quality"],
        turnover_arg=(
            f"--turnover_dir {os.path.join(OUTDIR, 'strain_turnover')}"
            if (GENE_DE_NOVO_PARAMS or {}).get("strain_gate") == "carrier"
            else ""
        ),
    shell:
        """
        alleleflux-gene-de-novo \
            --comparison {wildcards.timepoints}-{wildcards.groups} \
            --cache_dir {params.cache_dir} \
            --allele_analysis_dir {params.allele_analysis_dir} \
            --metadata {input.metadata} \
            --prodigal_fasta {input.prodigal_fasta} \
            --fasta {input.fasta} \
            --mag_mapping {input.mag_mapping} \
            --output_dir {params.output_dir} \
            {params.knob_args} \
            --min_cages {params.min_cages} \
            {params.fdr_arg} \
            --min_cov {params.min_cov} \
            --min_freq {params.min_freq} \
            --fdr {params.fdr} \
            --min_base_quality {params.min_base_quality} \
            {params.turnover_arg} \
            --cpus {threads}
        """
