import os
import pandas as pd
from glob import glob
from collections import defaultdict
from multiprocessing import Pool
import logging


# Load the configuration file
configfile: os.path.join(workflow.basedir, "config.yml")


# Define the global data_type variable to be used across all Snakemake files
DATA_TYPE = config["analysis"]["data_type"]

OUTDIR = config["output"]["root_dir"]

if DATA_TYPE == "single":
    OUTDIR = os.path.join(OUTDIR, "single_timepoint")
elif DATA_TYPE == "longitudinal":
    OUTDIR = os.path.join(OUTDIR, "longitudinal")

timepoints_labels = [
    tp[0] if len(tp) == 1 else f"{tp[0]}_{tp[1]}"
    for tp in config["analysis"]["timepoints_combinations"]
]  # ["T1", "T1_T2", "T2_T3"]

groups_labels = [
    f"{gr[0]}_{gr[1]}" for gr in config["analysis"]["groups_combinations"]
]  # ["G1_G2", "G2_G4"]

# Flatten the list of group values from the groups_combinations config.
group_values = sorted({item for gr in config["analysis"]["groups_combinations"] for item in gr})
# Build the regex: optionally an underscore and one of the allowed group values.
group_str_regex = "(_({}))?".format("|".join(group_values))


wildcard_constraints:
    groups=f"({'|'.join(groups_labels)})",
    timepoints=f"({'|'.join(timepoints_labels)})",
    taxon="(domain|phylum|class|order|family|genus|species)",
    test_type="(two_sample_unpaired|two_sample_paired|single_sample|lmm)",
    group_str=group_str_regex,


def get_mags_by_eligibility(timepoints, groups, eligibility_type):
    """
    Read the eligibility file for a given timepoint-group combination and return a list of MAG IDs.

    Parameters:
timepoints (str): The timepoints label
        groups (str): The groups label
        eligibility_type (str or None):
            - "two_sample_unpaired": only return MAG IDs where unpaired_test_eligible is True.
            - "two_sample_paired": only return MAG IDs where paired_test_eligible is True.
            - "lmm": only return MAG IDs where unpaired_test_eligible is True.
            - "all": return MAG IDs that are eligible for any of the tests.
    """
    eligibility_file = os.path.join(
        OUTDIR, f"eligibility_table_{timepoints}-{groups}.tsv"
    )
    df = pd.read_csv(eligibility_file, sep="\t")

    if eligibility_type == "two_sample_unpaired":
        return df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
    elif eligibility_type == "two_sample_paired":
        return df.loc[df["paired_test_eligible"] == True, "MAG_ID"].tolist()
    elif eligibility_type == "lmm":
        # LMM uses unpaired data, so we check unpaired eligibility.
        return df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
    # Return MAGs from all eligible columns
    elif eligibility_type == "all":
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
    else:
        raise ValueError(
            f"Unknown eligibility type: {eligibility_type}. "
            "Please use 'two_sample_unpaired', 'two_sample_paired', 'lmm' or 'all'."
        )



def get_single_sample_entries(timepoints, groups):
    """
    Reads the eligibility file and returns a list of tuples (MAG_ID, group)
    for each column matching 'single_sample_eligible_*' that evaluates to True.
    This allows a MAG to be eligible for multiple single-sample tests.
    """
    eligibility_file = os.path.join(
        OUTDIR, f"eligibility_table_{timepoints}-{groups}.tsv"
    )
    df = pd.read_csv(eligibility_file, sep="\t")
    sample_entries = []
    # Identify all columns with the generic pattern.
    sample_cols = [
        col for col in df.columns if col.startswith("single_sample_eligible_")
    ]

    # If there are no single_sample_eligible columns or all are NaN
    # (which happens for single data_type), return empty list
    if not sample_cols or df[sample_cols].isna().all().all():
        return []

    for _, row in df.iterrows():
        mag = row["MAG_ID"]
        for col in sample_cols:
            if pd.notna(row[col]) and row[col] == True:
                # Extract the sample group (e.g. "control", "fat") from the column name.
                group = col.replace("single_sample_eligible_", "")
                sample_entries.append((mag, group))
    return sample_entries
