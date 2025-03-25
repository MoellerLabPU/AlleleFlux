rule gene_scores:
    input:
        pvalue_table=os.path.join(
            OUTDIR,
            "significance_tests",
            "{test_type}_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}.tsv.gz",
        ),
    output:
        combined=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "gene_scores_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}_gene_scores_combined.tsv",
        ),
        individual=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "gene_scores_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}_gene_scores_individual.tsv",
        ),
        overlapping=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "gene_scores_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}_gene_scores_overlapping.tsv",
        ),
    params:
        prefix="{mag}_{test_type}{group_str}",
        pValue_threshold=config.get("p_value_threshold", 0.05),
        outDir=os.path.join(
            OUTDIR, "scores", "processed", "gene_scores_{timepoints}-{groups}"
        ),
        lmm_format=lambda wildcards: "--lmm_format" if wildcards.test_type == "lmm" else "",
    resources:
        time=config["time"]["general"],
    shell:
        """
        alleleflux-gene-scores \
            --pValue_table {input.pvalue_table} \
            --pValue_threshold {params.pValue_threshold} \
            --output_dir {params.outDir} \
            --prefix {params.prefix} \
            {params.lmm_format}
        """


rule detect_outlier_genes:
    input:
        mag_score=os.path.join(
            OUTDIR,
            "scores",
            "intermediate",
            "MAG_scores_{timepoints}-{groups}",
            "{mag}_score_{test_type}{group_str}.tsv",
        ),
        gene_scores=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "gene_scores_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}_gene_scores_individual.tsv",
        ),
    output:
        os.path.join(
            OUTDIR,
            "outlier_genes",
            "{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}_outlier_genes.tsv",
        ),
    params:
        # No params needed with entry points
    resources:
        time=config["time"]["general"],
    shell:
        """
        alleleflux-outliers \
            --mag_file {input.mag_score} \
            --mag_id {wildcards.mag} \
            --gene_file {input.gene_scores} \
            --out_fPath {output}
        """
