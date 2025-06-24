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

timepoints_labels = []
focus_timepoints = {}

for time_combo in config["analysis"]["timepoints_combinations"]:
    # Standardized format: all entries are dictionaries with "timepoint" key
    timepoint = time_combo["timepoint"]
    
    if len(timepoint) == 1 and DATA_TYPE == "single":
        # Single timepoint
        tp = timepoint[0]
        timepoints_labels.append(tp)
        # For single data type, the timepoint itself is the focus
        focus_timepoints[tp] = tp
    elif len(timepoint) == 2 and DATA_TYPE == "longitudinal":
        # Multiple timepoints (for longitudinal analysis)
        label = f"{timepoint[0]}_{timepoint[1]}"
        timepoints_labels.append(label)
        if "focus" in time_combo:
            # Ensure focus is one of the two timepoints
            if time_combo["focus"] in timepoint:
                focus_timepoints[label] = time_combo["focus"]
            else:
                raise ValueError(f"Invalid focus timepoint '{time_combo['focus']}' for combination '{label}'. Must be one of {timepoint}")
        else:
            # Default to second timepoint if focus not specified
            logging.warning(f"No focus specified for timepoint combination {label}. Using '{timepoint[1]}' as default.")
            focus_timepoints[label] = timepoint[1]

# Define valid focus timepoint values for each timepoint label
valid_focus_timepoints = {}
for tp in timepoints_labels:
    if "_" in tp and DATA_TYPE == "longitudinal":  # It's a two-timepoint combination with focus
        valid_focus_timepoints[tp] = tp.split("_")
    elif DATA_TYPE == "single":
        # For single timepoint, the timepoint itself is the focus
        valid_focus_timepoints[tp] = [tp]

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
    test_type="(two_sample_unpaired|two_sample_paired|single_sample|lmm|lmm_across_time|cmh|cmh_across_time|)",
    group_str=group_str_regex,
    # timepoint_str=timepoint_str_regex,
    # Constraint for focus timepoints - all possible values from the timepoint pairs
    focus_tp="|".join(set([tp for tp_pair in valid_focus_timepoints.values() for tp in tp_pair if tp])),
    
# Function to get sample information from metadata file
def get_sample_info():
    """
    Helper function to retrieve sample information.

    This function is responsible for loading and parsing sample metadata,
    typically from a configuration file or input source. It should return
    information necessary for processing samples in the workflow.

    Returns:
        dict or pandas.DataFrame: Sample information containing metadata
        required for the pipeline execution.
    """
    metadata_path = config["input"]["metadata_path"]
    
    if not metadata_path:
        raise ValueError("metadata_path must be provided in the config file")
    
    # Read metadata file
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    
    # Validate required columns
    required_cols = ["sample_id", "bam_path"]
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        raise ValueError(f"Metadata file is missing required columns: {', '.join(missing_cols)}")
    
    # Create mapping from sample ID to BAM path
    sample_to_bam_map = dict(zip(metadata_df["sample_id"], metadata_df["bam_path"]))
    sample_ids = list(metadata_df["sample_id"])
    
    # Validate that all BAM files exist
    for sample_id, bam_path in sample_to_bam_map.items():
        if not os.path.exists(bam_path):
            raise ValueError(f"BAM file not found for sample {sample_id}: {bam_path}")
    
    return sample_ids, sample_to_bam_map

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

    if eligibility_type == "two_sample_unpaired" or eligibility_type == "lmm":
        return df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
    elif eligibility_type == "two_sample_paired" or eligibility_type == "cmh":
        return df.loc[df["paired_test_eligible"] == True, "MAG_ID"].tolist()
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
            "Please use 'two_sample_unpaired', 'two_sample_paired', 'cmh', 'lmm' or 'all'."
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
