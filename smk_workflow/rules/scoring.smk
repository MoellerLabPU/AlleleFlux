rule significance_score_per_MAG:
    input:
        # {group_str} is empty for two-sample tests and "_{group}" for single-sample tests.
        pvalue_table=os.path.join(
            OUTDIR,
            "significance_tests",
            "{test_type}_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}.tsv.gz",
        ),
        gtdb_taxonomy=config["gtdb_file"],
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "intermediate",
            "MAG_scores_{timepoints}-{groups}",
            "{mag}_score_{test_type}{group_str}.tsv",
        ),
    params:
        group_by_column="MAG_ID",
        pValue_threshold=config.get("p_value_threshold", 0.05),
        lmm_format=lambda wildcards: "--lmm_format" if wildcards.test_type == "lmm" else "",
    resources:
        time=config["time"]["general"],
    shell:
        """
        alleleflux-scores \
            --gtdb_taxonomy {input.gtdb_taxonomy} \
            --pValue_table {input.pvalue_table} \
            --group_by_column {params.group_by_column} \
            --pValue_threshold {params.pValue_threshold} \
            --out_fPath {output} \
            {params.lmm_format}
        """

rule combine_MAG_scores:
    input:
        scores=lambda wc: expand(
            os.path.join(
                OUTDIR,
                "scores",
                    "intermediate",
                    "MAG_scores_{timepoints}-{groups}",
                    "{mag}_score_{test_type}{group_str}.tsv",
                ),
                timepoints=wc.timepoints,
                groups=wc.groups,
                test_type=wc.test_type,
                group_str=wc.group_str,
                mag=(
                    get_mags_by_eligibility(
                        wc.timepoints, wc.groups, eligibility_type=wc.test_type
                    )
                if wc.test_type != "single_sample"
                else [
                    mag
                    for mag, grp in get_single_sample_entries(wc.timepoints, wc.groups)
                    if f"_{grp}" == wc.group_str
                ]
            ),
        ),
    output:
        concatenated=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-MAGs.tsv",
        ),
    resources:
        time=config["time"]["general"],
    run:
        import logging
        import pandas as pd

        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
            level=logging.INFO,
        )

        dfs = []
        for file in input.scores:
            logging.info(f"Reading {file}")
            df = pd.read_csv(file, sep="\t")
            dfs.append(df)

        logging.info(
            f"Combining scores for {wildcards.timepoints}-{wildcards.groups} "
            f"({wildcards.test_type}{wildcards.group_str})"
        )

        combined_df = pd.concat(dfs, ignore_index=True)

        logging.info(f"Writing combined scores to {output.concatenated}")
        combined_df.to_csv(output.concatenated, sep="\t", index=False)

rule taxa_scores:
    input:
        concatenated=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-MAGs.tsv",
        ),
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-{taxon}.tsv",
        ),
    params:
        # No params needed with entry points
    resources:
        time=config["time"]["general"],
    shell:
        """
        alleleflux-taxa-scores \
            --input_df {input.concatenated} \
            --group_by_column {wildcards.taxon} \
            --out_fPath {output} 
        """
