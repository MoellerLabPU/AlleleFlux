import os
import pandas as pd
from glob import glob
from collections import defaultdict
from multiprocessing import Pool
import logging


# Load the configuration file
configfile: os.path.join(workflow.basedir, "config.yml")


OUTDIR = config["root_out"]

timepoints_labels = [
    f"{tp[0]}_{tp[1]}" for tp in config["timepoints_combinations"]
]  # ["T1_T2", "T2_T3"]
groups_labels = [
    f"{gr[0]}_{gr[1]}" for gr in config["groups_combinations"]
]  # ["G1_G2", "G2_G4"]

# Flatten the list of group values from the groups_combinations config.
group_values = sorted({item for gr in config["groups_combinations"] for item in gr})
# Build the regex: optionally an underscore and one of the allowed group values.
group_str_regex = "(_({}))?".format("|".join(group_values))


wildcard_constraints:
    groups=f"({'|'.join(groups_labels)})",
    timepoints=f"({'|'.join(timepoints_labels)})",
    taxon="(domain|phylum|class|order|family|genus|species)",
    test_type="(two_sample_unpaired|two_sample_paired|single_sample)",
    group_str=group_str_regex,


def get_mags_by_eligibility(timepoints, groups, eligibility_type=None):
    """
    Read the eligibility file for a given timepoint-group combination and return a list of MAG IDs.

    Parameters:
        eligibility_type (str or None):
            - "unpaired": only return MAG IDs where unpaired_test_eligible is True.
            - "paired": only return MAG IDs where paired_test_eligible is True.
            - None: return MAG IDs that are eligible by any of the criteria.
    """
    eligibility_file = os.path.join(OUTDIR, f"eligibility_table_{timepoints}-{groups}")
    df = pd.read_csv(eligibility_file, sep="\t")

    if eligibility_type == "unpaired":
        return df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
    elif eligibility_type == "paired":
        return df.loc[df["paired_test_eligible"] == True, "MAG_ID"].tolist()
    else:
        # Combine unpaired, paired, and any single-sample eligibility columns.
        single_cols = [
            col for col in df.columns if col.startswith("single_sample_eligible_")
        ]
        return (
            df[
                (df["unpaired_test_eligible"] == True)
                | (df["paired_test_eligible"] == True)
                | (df[single_cols].any(axis=1))
            ]["MAG_ID"]
            .unique()
            .tolist()
        )


def get_allele_analysis_targets():
    """
    Dynamically generate targets for allele analysis based on eligible MAGs.
    """
    targets = []
    for timepoints in timepoints_labels:
        for groups in groups_labels:
            # Get eligible MAGs for this timepoint-group combination
            eligible_mags = get_mags_by_eligibility(timepoints, groups)
            # Add targets for each eligible MAG
            for mag in eligible_mags:
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "allele_analysis",
                        f"allele_analysis_{timepoints}-{groups}",
                        f"{mag}_allele_frequency_changes_mean.tsv.gz",
                    )
                )
    return targets


# def get_two_sample_targets(test_type):
#     """
#     Generate targets for two-sample tests based on eligible MAGs.

#     Parameters:
#         test_type (str): either "unpaired" or "paired"
#     """
#     # Define subdirectory and file suffix based on the test type.
#     if test_type == "two_sample_unpaired":
#         subdir = "two_sample_unpaired"
#         suffix = "_two_sample_unpaired.tsv.gz"
#     elif test_type == "two_sample_paired":
#         subdir = "two_sample_paired"
#         suffix = "_two_sample_paired.tsv.gz"
#     else:
#         raise ValueError("test_type must be either 'unpaired' or 'paired'")

#     targets = []
#     for tp in timepoints_labels:
#         for gr in groups_labels:
#             for mag in get_mags_by_eligibility(tp, gr, eligibility_type=test_type):
#                 targets.append(
#                     os.path.join(
#                         OUTDIR,
#                         "significance_tests",
#                         f"{subdir}_{tp}-{gr}",
#                         f"{mag}{suffix}",
#                     )
#                 )
#     return targets


def get_single_sample_entries(timepoints, groups):
    """
    Reads the eligibility file and returns a list of tuples (MAG_ID, group)
    for each column matching 'single_sample_eligible_*' that evaluates to True.

    This allows a MAG to be eligible for multiple single-sample tests.
    """
    eligibility_file = os.path.join(OUTDIR, f"eligibility_table_{timepoints}-{groups}")
    df = pd.read_csv(eligibility_file, sep="\t")
    sample_entries = []
    # Identify all columns with the generic pattern.
    sample_cols = [
        col for col in df.columns if col.startswith("single_sample_eligible_")
    ]
    for _, row in df.iterrows():
        mag = row["MAG_ID"]
        for col in sample_cols:
            if row[col] == True:
                # Extract the sample group (e.g. "control", "fat") from the column name.
                group = col.replace("single_sample_eligible_", "")
                sample_entries.append((mag, group))
    return sample_entries


# def get_single_sample_targets():
#     """
#     Generate targets for single-sample tests. Each target is generated for a MAG and each
#     eligible sample group (e.g., control, fat) found in the eligibility file.
#     """
#     targets = []
#     for tp in timepoints_labels:
#         for gr in groups_labels:
#             for mag, group in get_single_sample_entries(tp, gr):
#                 targets.append(
#                     os.path.join(
#                         OUTDIR,
#                         "significance_tests",
#                         f"single_sample_{tp}-{gr}",
#                         f"{mag}_{group}_single_sample.tsv.gz",
#                     )
#                 )
#     return targets


def get_significance_scores_targets():
    targets = []
    # For two-sample tests (unpaired and paired)
    for tp in timepoints_labels:
        for gr in groups_labels:
            for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                for mag in get_mags_by_eligibility(tp, gr, eligibility_type=test_type):
                    targets.append(
                        os.path.join(
                            OUTDIR,
                            "scores",
                            "intermediate",
                            f"MAG_scores_{tp}-{gr}",
                            f"{mag}_score_{test_type}.tsv",
                        )
                    )
    # For single-sample tests
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag, group in get_single_sample_entries(tp, gr):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "scores",
                        "intermediate",
                        f"MAG_scores_{tp}-{gr}",
                        f"{mag}_score_single_sample_{group}.tsv",
                    )
                )
    return targets


def get_combined_scores_targets():
    targets = []
    # Two-sample tests: group_str is empty.
    for tp in timepoints_labels:
        for gr in groups_labels:
            for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                group_str = ""  # no group marker for two-sample tests
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "scores",
                        "processed",
                        "combined",
                        f"scores_{test_type}-{tp}-{gr}{group_str}-MAGs.tsv",
                    )
                )
    # Single-sample tests: group_str is "_" plus the sample group.
    for tp in timepoints_labels:
        for gr in groups_labels:
            # get_single_sample_entries returns (mag, group) pairs for a given timepoint and group.
            sample_entries = get_single_sample_entries(tp, gr)
            # print(sample_entries)
            unique_groups = sorted(set([grp for mag, grp in sample_entries]))
            # print(unique_groups)
            for grp in unique_groups:
                group_str = f"_{grp}"
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "scores",
                        "processed",
                        "combined",
                        f"scores_single_sample-{tp}-{gr}{group_str}-MAGs.tsv",
                    )
                )
    return targets


def get_taxa_scores_targets():
    targets = []
    tax_levels = ["phylum", "class", "order", "family", "genus", "species"]
    # For two-sample tests, group_str is empty.
    for tp in timepoints_labels:
        for gr in groups_labels:
            for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                group_str = ""  # no group marker for two-sample tests
                for taxon in tax_levels:
                    targets.append(
                        os.path.join(
                            OUTDIR,
                            "scores",
                            "processed",
                            "combined",
                            f"scores_{test_type}-{tp}-{gr}{group_str}-{taxon}.tsv",
                        )
                    )
    # For single-sample tests, group_str is "_" plus the sample group.
    for tp in timepoints_labels:
        for gr in groups_labels:
            sample_entries = get_single_sample_entries(tp, gr)
            unique_groups = sorted(set([grp for mag, grp in sample_entries]))
            for grp in unique_groups:
                group_str = f"_{grp}"
                for taxon in tax_levels:
                    targets.append(
                        os.path.join(
                            OUTDIR,
                            "scores",
                            "processed",
                            "combined",
                            f"scores_single_sample-{tp}-{gr}{group_str}-{taxon}.tsv",
                        )
                    )
    return targets


def get_gene_scores_targets():
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                group_str = ""  # no group marker for two-sample tests
                mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                for mag in mags:
                    prefix = f"{mag}_{test_type}{group_str}"
                    base_dir = os.path.join(
                        OUTDIR,
                        "scores",
                        "processed",
                        f"gene_scores_{tp}-{gr}",
                    )
                    targets.extend(
                        [
                            os.path.join(
                                base_dir, f"{prefix}_gene_scores_combined.tsv"
                            ),
                            os.path.join(
                                base_dir, f"{prefix}_gene_scores_individual.tsv"
                            ),
                            os.path.join(
                                base_dir, f"{prefix}_gene_scores_overlapping.tsv"
                            ),
                        ]
                    )

            sample_entries = get_single_sample_entries(tp, gr)
            for mag, grp in sample_entries:
                group_str = f"_{grp}"
                prefix = f"{mag}_single_sample{group_str}"
                base_dir = os.path.join(
                    OUTDIR,
                    "scores",
                    "processed",
                    f"gene_scores_{tp}-{gr}",
                )
                targets.extend(
                    [
                        os.path.join(base_dir, f"{prefix}_gene_scores_combined.tsv"),
                        os.path.join(base_dir, f"{prefix}_gene_scores_individual.tsv"),
                        os.path.join(base_dir, f"{prefix}_gene_scores_overlapping.tsv"),
                    ]
                )
    return targets


def get_outlier_gene_targets():
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                group_str = ""  # no group marker for two-sample tests
                mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                for mag in mags:
                    prefix = f"{mag}_{test_type}{group_str}"
                    base_dir = os.path.join(
                        OUTDIR,
                        "outlier_genes",
                        f"{tp}-{gr}",
                    )
                    targets.append(
                        os.path.join(base_dir, f"{prefix}_outlier_genes.tsv")
                    )

            sample_entries = get_single_sample_entries(tp, gr)
            for mag, grp in sample_entries:
                group_str = f"_{grp}"
                prefix = f"{mag}_single_sample{group_str}"
                base_dir = os.path.join(
                    OUTDIR,
                    "outlier_genes",
                    f"{tp}-{gr}",
                )
                targets.append(os.path.join(base_dir, f"{prefix}_outlier_genes.tsv"))
    return targets


localrules:
    all,


rule all:
    input:
        get_allele_analysis_targets(),
        # get_two_sample_targets("two_sample_unpaired"),
        # get_two_sample_targets("two_sample_paired"),
        # get_single_sample_targets(),
        # get_significance_scores_targets(),
        # get_combined_scores_targets(),
        get_taxa_scores_targets(),
        # get_gene_scores_targets(),
        get_outlier_gene_targets(),


rule analyze_alleles:
    input:
        mag_metadata_file=os.path.join(
            OUTDIR,
            "inputMetadata",
            "inputMetadata_{timepoints}-{groups}",
            "{mag}_metadata.tsv",
        ),
        eligibility_table=os.path.join(
            OUTDIR, "eligibility_table_{timepoints}-{groups}"
        ),
    output:
        # allele_changes=os.path.join(
        #     OUTDIR,
        #     "allele_analysis",
        #     "allele_analysis_{timepoints}-{groups}",
        #     "{mag}_allele_frequency_changes.tsv.gz",
        # ),
        mean_allele_changes=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean.tsv.gz",
        ),
    params:
        outDir=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
        ),
        scriptPath=config["scripts"]["analyze_alleles"],
        fasta=config["fasta"],
        breath_threshold=config.get("breath_threshold", 0.1),
        disable_zero_diff_filtering=(
            "--disable_zero_diff_filtering"
            if config.get("disable_zero_diff_filtering", False)
            else ""
        ),
    threads: config["cpus"]["analyze_alleles"]
    resources:
        mem_mb=config["memory"]["analyze_alleles"],
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --magID {wildcards.mag} \
            --mag_metadata_file {input.mag_metadata_file} \
            --fasta {params.fasta} \
            --breath_threshold {params.breath_threshold} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.disable_zero_diff_filtering}
            """


preprocess_enabled = config.get("preprocess_two_sample", False)
if preprocess_enabled:
    two_sample_input = os.path.join(
        OUTDIR,
        "significance_tests",
        "preprocessed_two_sample_{timepoints}-{groups}",
        "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz",
    )
else:
    two_sample_input = os.path.join(
        OUTDIR,
        "allele_analysis",
        "allele_analysis_{timepoints}-{groups}",
        "{mag}_allele_frequency_changes_mean.tsv.gz",
    )


rule preprocess_two_sample:
    input:
        mean_allele_changes=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean.tsv.gz",
        ),
    output:
        outPath=os.path.join(
            OUTDIR,
            "significance_tests",
            "preprocessed_two_sample_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean_preprocessed.tsv.gz",
        ),
    params:
        scriptPath=config["scripts"]["preprocess_two_sample"],
        alpha=config.get("alpha", 0.05),
        test_type=config.get("test_type", "t-test"),
    threads: config["cpus"]["quality_control"]
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --mean_changes_fPath {input.mean_allele_changes} \
            --cpus {threads} --alpha {params.alpha} \
            --output_fPath {output.outPath} \
            --test_type {params.test_type}
            """


rule two_sample_unpaired:
    input:
        # mean_allele_changes=os.path.join(
        #     OUTDIR,
        #     "allele_analysis",
        #     "allele_analysis_{timepoints}-{groups}",
        #     "{mag}_allele_frequency_changes_mean.tsv.gz",
        # ),
        mean_allele_changes=two_sample_input,
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_unpaired_{timepoints}-{groups}",
            "{mag}_two_sample_unpaired.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_unpaired_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["two_sample_unpaired"],
    threads: config["cpus"]["significance_test"]
    resources:
        # mem_mb=config["memory"]["significance_test"],
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --mean_changes_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir}
            """


rule two_sample_paired:
    input:
        # mean_allele_changes=os.path.join(
        #     OUTDIR,
        #     "allele_analysis",
        #     "allele_analysis_{timepoints}-{groups}",
        #     "{mag}_allele_frequency_changes_mean.tsv.gz",
        # ),
        mean_allele_changes=two_sample_input,
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "two_sample_paired_{timepoints}-{groups}",
            "{mag}_two_sample_paired.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "two_sample_paired_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["two_sample_paired"],
    threads: config["cpus"]["significance_test"]
    resources:
        # mem_mb=config["memory"]["significance_test"],
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --mean_changes_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --cpus {threads} \
            --output_dir {params.outDir}
            """


rule single_sample:
    input:
        mean_allele_changes=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean.tsv.gz",
        ),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "single_sample_{timepoints}-{groups}",
            "{mag}_single_sample_{group}.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "single_sample_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["single_sample"],
        max_zero_flag=(
            "--max_zero_count " + str(config["max_zero_count"])
            if config.get("max_zero_count", None) is not None
            else ""
        ),
    threads: config["cpus"]["significance_test"]
    resources:
        # mem_mb=config["memory"]["significance_test"],
        time=config["time"]["significance_test"],
    shell:
        """
        python {params.scriptPath} \
            --mean_changes_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --group {wildcards.group} \
            --cpus {threads} \
            --output_dir {params.outDir} \
            {params.max_zero_flag}
            """


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
        scriptPath=config["scripts"]["scores"],
        group_by_column="MAG_ID",
        pValue_threshold=config.get("p_value_threshold", 0.05),
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --gtdb_taxonomy {input.gtdb_taxonomy} \
            --pValue_table {input.pvalue_table} \
            --group_by_column {params.group_by_column} \
            --pValue_threshold {params.pValue_threshold} \
            --out_fPath {output}
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
        scriptPath=config["scripts"]["taxa_scores"],
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --input_df {input.concatenated} \
            --group_by_column {wildcards.taxon} \
            --out_fPath {output} 
        """


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
        scriptPath=config["scripts"]["gene_scores"],  # Update with the actual path to your script
        prefix="{mag}_{test_type}{group_str}",
        pValue_threshold=config.get("p_value_threshold", 0.05),
        outDir=os.path.join(
            OUTDIR, "scores", "processed", "gene_scores_{timepoints}-{groups}"
        ),
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --pValue_table {input.pvalue_table} \
            --pValue_threshold {params.pValue_threshold} \
            --output_dir {params.outDir} \
            --prefix {params.prefix} 
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
        scriptPath=config["scripts"]["outlier_detection"],
    resources:
        time=config["time"]["general"],
    shell:
        """
        python {params.scriptPath} \
            --mag_file {input.mag_score} \
            --mag_id {wildcards.mag} \
            --gene_file {input.gene_scores} \
            --out_fPath {output}
        """
