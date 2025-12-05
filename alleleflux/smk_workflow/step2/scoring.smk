

import pandas as pd
# import logging
# from alleleflux.scripts.utilities.logging_config import setup_logging
from snakemake.logging import logger

rule significance_score_per_MAG_standard:
    input:
        # Standard test types use a single pvalue table
        pvalue_table=os.path.join(
            OUTDIR,
            "significance_tests",
            "{test_type}_{timepoints}-{groups}",
            "{mag}_{test_type}{group_str}.tsv.gz",
        ),
        gtdb_taxonomy=config["input"]["gtdb_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "intermediate",
            "MAG_scores_{timepoints}-{groups}",
            "{mag}_score_{test_type}{group_str}.tsv"
        ),
    params:
        group_by_column="MAG_ID",
        pValue_threshold=config["statistics"].get("p_value_threshold", 0.05),
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-scores \
            --gtdb_taxonomy {input.gtdb_taxonomy} \
            --pValue_table {input.pvalue_table} \
            --group_by_column {params.group_by_column} \
            --pValue_threshold {params.pValue_threshold} \
            --out_fPath {output} \
            --mag_mapping_file {input.mag_mapping}
        """


rule significance_score_per_MAG_cmh:
    input:
        pvalue_table=os.path.join(
            OUTDIR,
            "significance_tests",
            "cmh_{timepoints}-{groups}",
            "{mag}_cmh.tsv.gz",
        ),
        gtdb_taxonomy=config["input"]["gtdb_path"],
        mag_mapping=config["input"]["mag_mapping_path"],
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "intermediate",
            "MAG_scores_{timepoints}-{groups}",
            "{mag}_score_cmh_{focus_tp}.tsv"
        ),
    params:
        # CMH-specific parameters
        tp1_name=lambda wildcards: wildcards.timepoints.split("_")[0] if DATA_TYPE == "longitudinal" else "",
        tp2_name=lambda wildcards: wildcards.timepoints.split("_")[1] if DATA_TYPE == "longitudinal" else "",
        pValue_threshold=config["statistics"].get("p_value_threshold", 0.05),
        group_by_column="MAG_ID",
        data_type=DATA_TYPE,
    resources:
        time=config["resources"]["time"]["general"],
    run:
        if params.data_type == "single":
            shell(
            """
            alleleflux-scores \
                --gtdb_taxonomy {input.gtdb_taxonomy} \
                --pValue_table {input.pvalue_table} \
                --group_by_column {params.group_by_column} \
                --pValue_threshold {params.pValue_threshold} \
                --out_fPath {output} \
                --mag_mapping_file {input.mag_mapping}
            """)
        elif params.data_type == "longitudinal":
            shell(
            """
            alleleflux-cmh-scores \
                --combined-file {input.pvalue_table} \
                --tp1-name {params.tp1_name} \
                --tp2-name {params.tp2_name} \
                --focus {wildcards.focus_tp} \
                --gtdb_taxonomy {input.gtdb_taxonomy} \
                --mag_id {wildcards.mag} \
                --threshold {params.pValue_threshold} \
                --out_fPath {output}
            """)



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
                if wc.test_type not in {"single_sample", "lmm_across_time", "cmh_across_time"}
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
            "MAG",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-MAGs.tsv",
        ),
    resources:
        time=config["resources"]["time"]["general"],
    run:
        # setup_logging()
        # logger = logging.getLogger(__name__)

        dfs = []
        for file in input.scores:
            logger.info(f"Reading {file}")
            df = pd.read_csv(file, sep="\t")
            dfs.append(df)

        logger.info(
            f"Combining scores for {wildcards.timepoints}-{wildcards.groups} "
            f"({wildcards.test_type}{wildcards.group_str})"
        )

        combined_df = pd.concat(dfs, ignore_index=True)

        logger.info(f"Writing combined scores to {output.concatenated}")
        combined_df.to_csv(output.concatenated, sep="\t", index=False)


rule combine_MAG_scores_cmh:
    input:
        scores=lambda wc: expand(
            os.path.join(
                OUTDIR,
                "scores",
                "intermediate",
                "MAG_scores_{timepoints}-{groups}",
                "{mag}_score_cmh_{focus_tp}.tsv",
            ),
            timepoints=wc.timepoints,
            groups=wc.groups,
            mag=get_mags_by_eligibility(wc.timepoints, wc.groups, eligibility_type="cmh"),
            focus_tp=[wc.focus_tp],
        ),
    output:
        concatenated=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "MAG",
            "scores_cmh-{timepoints}-{groups}-MAGs-{focus_tp}.tsv",
        ),
    resources:
        time=config["resources"]["time"]["general"],
    run:
        # setup_logging()
        # logger = logging.getLogger(__name__)

        dfs = []
        for file in input.scores:
            logger.info(f"Reading {file}")
            df = pd.read_csv(file, sep="\t")
            # Verify that the focus timepoint in the file matches the expected one
            if "focus_timepoint" in df.columns and df["focus_timepoint"].iloc[0] != wildcards.focus_tp:
                raise ValueError(f"Mismatched focus timepoint in {file}: "
                                 f"expected {wildcards.focus_tp}, found {df['focus_timepoint'].iloc[0]}")
            dfs.append(df)

        logger.info(
            f"Combining CMH scores for {wildcards.timepoints}-{wildcards.groups} "
            f"with focus timepoint {wildcards.focus_tp}"
        )

        if not dfs:
            raise ValueError(f"No valid CMH score files found for focus timepoint {wildcards.focus_tp}")

        combined_df = pd.concat(dfs, ignore_index=True)

        logger.info(f"Writing combined scores to {output.concatenated}")
        combined_df.to_csv(output.concatenated, sep="\t", index=False)

rule taxa_scores:
    input:
        concatenated=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "MAG",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-MAGs.tsv",
        ),
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "{taxon}",
            "scores_{test_type}-{timepoints}-{groups}{group_str}-{taxon}.tsv",
        ),
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-taxa-scores \
            --input_df {input.concatenated} \
            --group_by_column {wildcards.taxon} \
            --out_fPath {output} 
        """


rule taxa_scores_cmh:
    input:
        concatenated=os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "MAG",
            "scores_cmh-{timepoints}-{groups}-MAGs-{focus_tp}.tsv",
        ),
    output:
        os.path.join(
            OUTDIR,
            "scores",
            "processed",
            "combined",
            "{taxon}",
            "scores_cmh-{timepoints}-{groups}-{taxon}-{focus_tp}.tsv",
        ),
    resources:
        time=config["resources"]["time"]["general"],
    shell:
        """
        alleleflux-taxa-scores \
            --input_df {input.concatenated} \
            --group_by_column {wildcards.taxon} \
            --out_fPath {output} 
        """
