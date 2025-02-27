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


wildcard_constraints:
    groups="[a-zA-Z0-9_-]+",  # Allow alphanumeric, underscores, and hyphens
    timepoints="[a-zA-Z0-9_-]+",
    taxon="[a-zA-Z]+",
    test_type="[a-zA-Z_]+",


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


def get_two_sample_targets(test_type):
    """
    Generate targets for two-sample tests based on eligible MAGs.

    Parameters:
        test_type (str): either "unpaired" or "paired"
    """
    # Define subdirectory and file suffix based on the test type.
    if test_type == "unpaired":
        subdir = "two_sample_unpaired"
        suffix = "_two_sample_unpaired.tsv.gz"
    elif test_type == "paired":
        subdir = "two_sample_paired"
        suffix = "_two_sample_paired.tsv.gz"
    else:
        raise ValueError("test_type must be either 'unpaired' or 'paired'")

    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag in get_mags_by_eligibility(tp, gr, eligibility_type=test_type):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "significance_tests",
                        f"{subdir}_{tp}-{gr}",
                        f"{mag}{suffix}",
                    )
                )
    return targets


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


def get_single_sample_targets():
    """
    Generate targets for single-sample tests. Each target is generated for a MAG and each
    eligible sample group (e.g., control, fat) found in the eligibility file.
    """
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag, group in get_single_sample_entries(tp, gr):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "significance_tests",
                        f"single_sample_{tp}-{gr}",
                        f"{mag}_{group}_single_sample.tsv.gz",
                    )
                )
    return targets


localrules:
    all,


rule all:
    input:
        get_allele_analysis_targets(),
        get_two_sample_targets("unpaired"),
        get_two_sample_targets("paired"),
        get_single_sample_targets(),


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


rule two_sample_unpaired:
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
            "{mag}_{group}_single_sample.tsv.gz",
        ),
    params:
        min_sample_num=config["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "single_sample_{timepoints}-{groups}"
        ),
        scriptPath=config["scripts"]["single_sample"],
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
            --output_dir {params.outDir}
            """
