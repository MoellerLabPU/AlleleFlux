"""Regional contrast analysis rule.

This module runs per-MAG regional contrast analysis, which detects genes or
sliding windows with differential allele-frequency evolution between treatment
and control groups across paired hosts.

The analysis operates directly on the ``*_allele_frequency_changes_mean.tsv.gz``
file produced by allele analysis (no statistical preprocessing required), and is
therefore only applicable to longitudinal data.

Rules are independent per MAG and run in parallel across the DAG.
"""

import os


rule regional_contrast:
    """Detect genomic regions with differential allele-frequency evolution.

    For each MAG, computes Total-Variation distance site scores, aggregates
    them into gene and/or sliding-window regions per host × group, ranks
    scores into genome-wide percentiles, contrasts treatment vs control, and
    tests the contrast across hosts with a Wilcoxon signed-rank test and a
    one-sample t-test (both directional), followed by Benjamini–Hochberg FDR
    correction.

    Runs independently per MAG; Snakemake will parallelise across MAGs.
    """
    input:
        get_allele_analysis_input_path(gr_wildcard="{treatment}_{control}"),
    output:
        per_host=os.path.join(
            OUTDIR,
            "regional_contrast",
            "regional_contrast_{timepoints}-{treatment}_{control}",
            "{mag}_regional_contrast_per_host_region.tsv.gz",
        ),
        summary=os.path.join(
            OUTDIR,
            "regional_contrast",
            "regional_contrast_{timepoints}-{treatment}_{control}",
            "{mag}_regional_contrast_region_summary.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTDIR,
            "regional_contrast",
            "regional_contrast_{timepoints}-{treatment}_{control}",
        ),
        prefix="{mag}_regional_contrast",
        mode=config.get("regional_contrast", {}).get("mode", "both"),
        window_size=config.get("regional_contrast", {}).get("window_size", 1000),
        agg_method=config.get("regional_contrast", {}).get("agg_method", "median"),
        trim_fraction=(
            "--trim_fraction "
            + str(config.get("regional_contrast", {}).get("trim_fraction", 0.1))
            if config.get("regional_contrast", {}).get("agg_method", "median")
            == "trimmed_mean"
            else ""
        ),
        min_informative_sites=config.get("regional_contrast", {}).get(
            "min_informative_sites", 5
        ),
        min_informative_fraction=config.get("regional_contrast", {}).get(
            "min_informative_fraction", 0.0
        ),
        use_fisher=(
            "--use_fisher"
            if config.get("regional_contrast", {}).get("use_fisher", True)
            else ""
        ),
    resources:
        mem_mb=get_mem_mb("regional_contrast"),
        time=get_time("regional_contrast"),
    shell:
        """
        alleleflux-regional-contrast \
            --input {input} \
            --output_dir {params.output_dir} \
            --prefix {params.prefix} \
            --treatment_group {wildcards.treatment} \
            --control_group {wildcards.control} \
            --mode {params.mode} \
            --window_size {params.window_size} \
            --agg_method {params.agg_method} \
            {params.trim_fraction} \
            --min_informative_sites {params.min_informative_sites} \
            --min_informative_fraction {params.min_informative_fraction} \
            {params.use_fisher}
        """
