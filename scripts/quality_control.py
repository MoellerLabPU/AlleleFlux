#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from glob import glob
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from utilities import calculate_mag_sizes, load_mag_metadata_file


def init_worker(metadata, mag_sizes):
    """
    Initialize worker process with metadata and MAG sizes.

    This function sets up global dictionaries for metadata and MAG sizes
    that can be accessed by worker processes.

    Parameters:
        metadata (dict): A dictionary containing metadata information.
        mag_sizes (dict): A dictionary containing sizes of MAGs (Metagenome-Assembled Genomes).

    Returns:
        None
    """
    global metadata_dict
    global mag_size_dict
    metadata_dict = metadata
    mag_size_dict = mag_sizes


def process_mag_files(args):
    sample_id, profile_fPath, mag_id, breadth_threshold = args

    # Build a dictionary to store coverage stats
    result = {
        "sample_id": sample_id,
        "MAG_ID": mag_id,
        "file_path": profile_fPath,
        "group": metadata_dict[sample_id]["group"],
        "subjectID": metadata_dict[sample_id]["subjectID"],
        "time": metadata_dict[sample_id]["time"],
        "replicate": metadata_dict[sample_id]["replicate"],
        "genome_size": mag_size_dict.get(mag_id, 0),
        "breadth": None,
        "breadth_threshold_passed": False,
        "breadth_fail_reason": "",
    }

    # Check if MAG size is known
    mag_size = mag_size_dict.get(mag_id)
    if mag_size is None or mag_size <= 0:
        msg = f"MAG size for {mag_id} not found or invalid."
        logging.error(f"{msg} sample={sample_id}, profile_fPath={profile_fPath}")
        result["breadth_fail_reason"] = msg
        return result

    df = pd.read_csv(profile_fPath, sep="\t", dtype={"gene_id": str})

    # Compute breadth coverage
    positions_with_coverage = df[df["total_coverage"] >= 1].shape[0]
    breadth = positions_with_coverage / mag_size
    result["breadth"] = breadth

    if breadth < breadth_threshold:
        msg = f"Breadth {breadth:.2%} < threshold {breadth_threshold:.2%}"
        logging.info(f"{msg} - sample={sample_id}, mag={mag_id}")
        result["breadth_fail_reason"] = msg
    else:
        result["breadth_threshold_passed"] = True
    return result


def check_timepoints(df):
    # Initialize columns
    df["two_timepoints_passed"] = False

    # Get samples that passed coverage threshold
    passed_samples = df[df["breadth_threshold_passed"]]

    # Check that there are exactly 2 unique timepoints globally.
    unique_times = passed_samples["time"].unique()
    if len(unique_times) > 2:
        raise ValueError(f"More than 2 unique timepoints found: {unique_times}")

    if passed_samples.empty:
        logging.warning("No samples passed coverage threshold")
        return df

    # Find subjects with MAG in exactly 2 timepoints
    valid_subjects = (
        passed_samples.groupby("subjectID")["time"]
        .nunique()
        .where(lambda x: x == 2)
        .dropna()
        .index
    )

    # Mark passed samples from valid subjects
    valid_mask = passed_samples["subjectID"].isin(valid_subjects)
    df.loc[passed_samples[valid_mask].index, "two_timepoints_passed"] = True

    # Add debug info
    if not valid_subjects.empty:
        logging.info(f"Subjects with 2 valid timepoints: {len(valid_subjects)}")
    else:
        logging.warning("No subjects have the MAG in 2 timepoints")

    return df


def add_subject_count_per_group(df):
    # Get samples that passed timepoint validation
    valid_subjects = df[df["two_timepoints_passed"]]

    # Check that there are exactly 2 unique groups in the valid samples.
    unique_groups = valid_subjects["group"].unique()
    if len(unique_groups) > 2:
        raise ValueError(f"More than 2 groups found: {unique_groups}")

    if valid_subjects.empty:
        df["subjects_per_group"] = np.nan
        df["replicates_per_group"] = np.nan
        return df

    # Count unique subjects per group in valid samples
    group_counts = (
        valid_subjects.groupby("group")["subjectID"]
        .nunique()  # Count unique subjects per group
        .reset_index(name="subjects_per_group")
    )

    # Count replicates per group in valid samples
    replicate_counts = (
        valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="replicates_per_group")
    )

    # Merge both counts into one counts DataFrame keyed by group.
    counts = group_counts.merge(replicate_counts, on="group", how="outer")

    # Merge the counts back into the original DataFrame based on group
    df = df.merge(counts, on="group", how="left")

    # For rows where two_timepoints_passed is False, set the counts to NaN.
    df.loc[
        ~df["two_timepoints_passed"], ["subjects_per_group", "replicates_per_group"]
    ] = np.nan

    # Log summary
    logging.info("Subject and replicate counts per group (valid samples):")
    logging.info(f"\n{counts}")
    return df


def count_paired_replicates(df):
    valid_subjects = df[df["two_timepoints_passed"]]

    # Check that there are exactly 2 unique groups in the valid samples.
    unique_groups = valid_subjects["group"].unique()
    if len(unique_groups) > 2:
        raise ValueError(f"More than 2 groups found: {unique_groups}")

    if valid_subjects.empty:
        df["subjects_per_group"] = np.nan
        df["replicates_per_group"] = np.nan
        return df

    rep_group_count = (
        valid_subjects.groupby("replicate")["group"]
        .nunique()
        .reset_index(name="group_count")
    )
    paired_reps = rep_group_count[rep_group_count["group_count"] == 2]["replicate"]
    paired_valid_subjects = valid_subjects[
        valid_subjects["replicate"].isin(paired_reps)
    ]

    # Count the number of unique paired replicates per group
    paired_counts = (
        paired_valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="paired_replicates_per_group")
    )

    # Merge the paired replicate counts back into the original DataFrame based on 'group'
    df = df.merge(paired_counts, on="group", how="left")

    # Set counts to NaN where:
    # - Samples failed timepoint check OR
    # - Samples passed timepoint check but replicate isn't paired
    valid_mask = df["two_timepoints_passed"]
    paired_mask = df["replicate"].isin(paired_reps)
    df["paired_replicates_per_group"] = np.where(
        valid_mask & paired_mask, df["paired_replicates_per_group"], np.nan
    )

    # Log summary
    logging.info("Replicates paired per group")
    logging.info(f"\n{paired_counts}")
    return df


def process_mag(args):
    mag_id, metadata_file, mag_size_dict, output_dir, breadth_threshold, cpus = args
    logging.info(f"Processing MAG: {mag_id}")

    metadata_dict, sample_files_with_mag_id = load_mag_metadata_file(
        metadata_file, mag_id, breadth_threshold
    )

    if not sample_files_with_mag_id:
        raise ValueError(
            f"No samples found in the metadata file {metadata_file}. For MAG: {mag_id}. Exiting."
        )

    num_procs = min(cpus, len(sample_files_with_mag_id))
    logging.info(
        f"Processing {len(sample_files_with_mag_id)} samples for {mag_id} with {num_procs} processes."
    )

    # Load sample data in parallel
    with Pool(
        processes=num_procs,
        initializer=init_worker,
        initargs=(metadata_dict, mag_size_dict),
    ) as pool:
        results_list = list(
            pool.imap_unordered(process_mag_files, sample_files_with_mag_id),
        )
    # Build results DataFrame
    df_results = pd.DataFrame(results_list)

    # Check that each subject has exactly 2 timepoints.
    df_results = check_timepoints(df_results)

    # Count the number of unique subjects and replicates per group
    df_results = add_subject_count_per_group(df_results)

    # Count the number of paired replicates per group
    df_results = count_paired_replicates(df_results)

    out_file = os.path.join(output_dir, f"{mag_id}_QC.tsv")
    df_results.to_csv(out_file, sep="\t", index=False)
    logging.info(f"Saved QC report for {mag_id} to {out_file}")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Script to detect which samples for a MAG pass the breadth threshold and output a pass/fail table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # parser.add_argument("--magID", required=True, help="MAG ID to process")
    parser.add_argument(
        "--rootDir",
        required=True,
        help="Path to the root directory containing metadata files.",
    )
    parser.add_argument(
        "--fasta", required=True, help="FASTA file with contigs (for MAG size)"
    )
    parser.add_argument(
        "--breadth_threshold",
        type=float,
        default=0.1,
        help="Minimum breadth coverage required to pass.",
    )
    parser.add_argument(
        "--cpus", type=int, default=cpu_count(), help="Number of processors to use."
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to write output table."
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    mag_size_dict = calculate_mag_sizes(args.fasta)

    metadata_files = glob(os.path.join(args.rootDir, "*_metadata.tsv"))
    if not metadata_files:
        logging.error("No *_metadata.tsv files found in input directory")
        sys.exit(1)

    tasks = [
        (
            os.path.basename(metadata_file).split("_metadata")[
                0
            ],  # mag_id extracted from filename
            metadata_file,
            mag_size_dict,
            args.output_dir,
            args.breadth_threshold,
            args.cpus,
        )
        for metadata_file in metadata_files
    ]

    for task in tasks:
        process_mag(task)

    logging.info("QC analysis completed.")


if __name__ == "__main__":
    main()
