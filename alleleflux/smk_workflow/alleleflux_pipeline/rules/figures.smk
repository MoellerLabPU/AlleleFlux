"""Score figures per comparison and the one-page HTML run report.

``plot_scores`` runs once per {timepoints}-{groups} comparison on that
comparison's combined MAG score tables plus the run-once significant-sites
summary (the smallest q per MAG and cell), and writes the heatmap and the two
boxplots per test family under ``figures/<comparison>/`` with a done-file
listing them.  ``run_report`` then runs once over every comparison's done-file,
the eligibility tables and the site-depth-filter summary, and writes
``figures/run_report.html`` with every figure embedded.

Both are ``localrules``: seconds of matplotlib, no profiles read, so a SLURM
submission would cost more than the job.  Targets come from
``generate_figure_targets()`` (dynamic_targets.smk), which reads the
eligibility checkpoints through generate_mag_scores_targets, so the rules'
inputs are lambdas (evaluated during DAG resolution, not at parse time).
"""

import os


localrules:
    plot_scores,
    run_report,


def _cell_stats_path():
    return os.path.join(
        OUTDIR, "p_value_summary", "significant_sites_summary",
        "significant_sites_mag_cell_stats_long.tsv",
    )


rule plot_scores:
    input:
        scores=lambda w: generate_mag_scores_targets(w.timepoints, w.groups),
        cell_stats=_cell_stats_path(),
        # the tree, when configured, is an input so a new tree redraws the figures
        tree=[FIGURES["tree_path"]] if FIGURES["tree_path"] else [],
    output:
        done=os.path.join(OUTDIR, "figures", "{timepoints}-{groups}", "score_figures.done"),
    params:
        output_dir=os.path.join(OUTDIR, "figures", "{timepoints}-{groups}"),
        # the two group labels in config order (a label may contain underscores)
        groups=lambda w: " ".join(groups_for_label(w.groups)),
        families_arg=(
            "--families " + " ".join(FIGURES["families"]) if FIGURES["families"] else ""
        ),
        tree_arg=f"--tree \"{FIGURES['tree_path']}\"" if FIGURES["tree_path"] else "",
        fmt=FIGURES["format"],
        min_mags_per_phylum=FIGURES["min_mags_per_phylum"],
        # stars at the same threshold the pipeline calls a site significant
        q_threshold=config["statistics"].get("p_value_threshold", 0.05),
        data_type=DATA_TYPE,
        run_name=config.get("run_name", ""),
    shell:
        """
        alleleflux-plot-scores \
            --run_dir "{OUTDIR}" \
            --comparison {wildcards.timepoints}-{wildcards.groups} \
            --groups {params.groups} \
            --output_dir "{params.output_dir}" \
            --format {params.fmt} \
            --q_threshold {params.q_threshold} \
            --min_mags_per_phylum {params.min_mags_per_phylum} \
            --data_type {params.data_type} \
            --run_name "{params.run_name}" \
            {params.families_arg} {params.tree_arg}
        """


def _figure_sentinels(wildcards):
    """Every comparison's done-file (the report embeds all of their figures)."""
    return [t for t in generate_figure_targets() if t.endswith("score_figures.done")]


def _eligibility_tables(wildcards):
    return [
        os.path.join(OUTDIR, f"eligibility_table_{tp}-{gr}.tsv")
        for tp in timepoints_labels
        for gr in groups_labels
    ]


rule run_report:
    input:
        sentinels=_figure_sentinels,
        eligibility=_eligibility_tables,
        cell_stats=_cell_stats_path(),
        site_filter=lambda w: generate_site_depth_filter_summary_targets(),
    output:
        report=os.path.join(OUTDIR, "figures", "run_report.html"),
    params:
        figures_dir=os.path.join(OUTDIR, "figures"),
        comparisons=" ".join(f"{tp}-{gr}" for tp in timepoints_labels for gr in groups_labels),
        # the per-MAG metadata index lives under the reuse dir when profiles were reused
        metadata_dir=os.path.join(REUSE_DIR, "inputMetadata"),
        metadata_path=config["input"]["metadata_path"],
        q_threshold=config["statistics"].get("p_value_threshold", 0.05),
        data_type=DATA_TYPE,
        permuted_arg="--permuted" if (GENERATE_PERMUTED or PERMUTED_METADATA_DIR) else "",
    shell:
        """
        alleleflux-run-report \
            --run_dir "{OUTDIR}" \
            --figures_dir "{params.figures_dir}" \
            --output "{output.report}" \
            --comparisons {params.comparisons} \
            --metadata_dir "{params.metadata_dir}" \
            --metadata_path "{params.metadata_path}" \
            --q_threshold {params.q_threshold} \
            --data_type {params.data_type} \
            {params.permuted_arg}
        """
