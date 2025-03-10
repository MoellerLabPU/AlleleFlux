import argparse
import gc
import logging
import os
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import stats
from tqdm import tqdm
from utilities import calculate_mag_sizes, load_mag_metadata_file

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]


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
    """
    Processes a MAG (Metagenome-Assembled Genome) file and adds metadata information.

    Parameters:
        args (tuple): A tuple containing the following elements:
            sample_id (str): The sample identifier.
            filepath (str): The path to the MAG file.
            mag_id (str): The MAG identifier.
            breath_threshold (float): The threshold for breadth coverage.

    Returns:
        pd.DataFrame or None: A DataFrame with added metadata and calculated breadth if the breadth is above the threshold,
                              otherwise None if the breadth is below the threshold or if the MAG size is not found.

    Notes:
        - The function reads the MAG file from the given filepath.
        - Adds sample_id, group, and subjectID columns to the DataFrame.
        - Inserts the MAG_ID as the first column.
        - Calculates the breadth of coverage for the MAG.
        - If the breadth is below the given threshold, the function logs a message and returns None.
        - If the breadth is above the threshold, the function adds the breadth and genome size to the DataFrame and returns it.
    """
    sample_id, filepath, mag_id, breath_threshold = args
    df = pd.read_csv(filepath, sep="\t", dtype={"gene_id": str})

    # Add sample_id column
    df["sample_id"] = sample_id
    # Add metadata columns
    metadata_info = metadata_dict.get(sample_id)
    df["group"] = metadata_info["group"]
    df["subjectID"] = metadata_info["subjectID"]
    df["time"] = metadata_info["time"]
    df["replicate"] = metadata_info["replicate"]
    # This adds MAG_ID as the first column
    df.insert(0, "MAG_ID", mag_id)

    # Get MAG size (total number of positions in the MAG)
    mag_size = mag_size_dict.get(mag_id)
    if mag_size is None:
        logging.warning(f"Size for MAG {mag_id} not found in sample {sample_id}.")
        return None  # Skip this sample-MAG combination

    # Calculate the number of positions with total_coverage >= 1
    positions_with_coverage = df[df["total_coverage"] >= 1].shape[0]

    # Calculate breadth
    breadth = positions_with_coverage / mag_size

    if breadth < breath_threshold:
        logging.info(
            f"MAG {mag_id} in sample {sample_id} has breadth {breadth:.2%}, which is less than {breath_threshold:.2%}. Skipping this sample-MAG combination."
        )
        return None  # Skip this sample-MAG combination
    else:
        df["breadth"] = breadth
        df["genome_size"] = mag_size
        df = calculate_frequencies(df, mag_id)
        return df


def calculate_frequencies(df, mag_id):
    """
    Calculate nucleotide frequencies for a given DataFrame.

    This function takes a DataFrame containing nucleotide counts and total coverage,
    and calculates the frequency of each nucleotide (A, C, T, G) by dividing the count
    of each nucleotide by the total coverage. The resulting frequencies are added as new
    columns to the DataFrame.

    Parameters:
        df (pandas.DataFrame): A DataFrame with columns "A", "C", "T", "G", and "total_coverage".

    Returns:
        pandas.DataFrame: The input DataFrame with additional columns for nucleotide frequencies:
                          "A_frequency", "C_frequency", "T_frequency", and "G_frequency".
    """
    total_coverage = df["total_coverage"]
    # logging.info(f"Calculating nucleotide frequencies for MAG {mag_id}.")
    # Calculate frequencies directly
    df["A_frequency"] = df["A"] / total_coverage
    df["T_frequency"] = df["T"] / total_coverage
    df["G_frequency"] = df["G"] / total_coverage
    df["C_frequency"] = df["C"] / total_coverage

    return df


def save_nucleotide_frequencies(data_dict, output_dir, mag_id):
    mag_df = pd.concat(
        [df for subject_dict in data_dict.values() for df in subject_dict.values()],
        ignore_index=True,
    )

    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Saving nucleotide frequencies for MAG {mag_id} to {output_dir}")
    mag_df.to_csv(
        os.path.join(output_dir, f"{mag_id}_nucleotide_frequencies.tsv.gz"),
        index=False,
        sep="\t",
        compression="gzip",
    )


def create_data_dict(data_list):
    """
    Organizes a list of DataFrames into a nested dictionary structure based on subjectID and timepoint.

    Parameters:
        data_list (list of pandas.DataFrame): A list of DataFrames, each containing columns 'subjectID', 'time', and 'sample_id'.

    Returns:
        dict: A nested dictionary where the first level keys are subjectIDs and the second level keys are timepoints.
              The values are the corresponding DataFrames.

    Raises:
        ValueError: If a DataFrame contains multiple unique subjectIDs or timepoints.
    """
    # Organize DataFrames by subjectID
    data_dict = {}
    for df in data_list:
        # Validate and extract subjectID
        subject_ids = df["subjectID"].unique()
        if len(subject_ids) != 1:
            raise ValueError(
                f"Multiple subjectIDs found in DataFrame for sample {df['sample_id'].iloc[0]}"
            )
        subjectID = subject_ids[0]
        # Validate and extract group
        timepoints = df["time"].unique()
        if len(timepoints) != 1:
            raise ValueError(
                f"Multiple timepoints found in DataFrame for sample {df['sample_id'].iloc[0]}"
            )
        timepoint = timepoints[0]

        # Store the DataFrame for each subjectID and timepoint
        if subjectID not in data_dict:
            data_dict[subjectID] = {}
        data_dict[subjectID][timepoint] = df
    return data_dict


def calculate_allele_frequency_changes(data_dict, output_dir, mag_id):
    logging.info("Identifying unique timepoints.")

    unique_timepoints = set()
    for subject_data in data_dict.values():
        unique_timepoints.update(subject_data.keys())
    if len(unique_timepoints) != 2:
        raise ValueError(
            f"Expected exactly 2 unique timepoints, found {len(unique_timepoints)}."
        )

    # Unpack the two timepoints
    timepoint_1, timepoint_2 = unique_timepoints

    logging.info(
        f"Calculating change in allele frequency between {timepoint_1} and {timepoint_2} for each position between the same subjectID."
    )
    # Get sets of subjectIDs in each timepoint
    subjectIDs_timepoint1 = {
        subjectID for subjectID in data_dict if timepoint_1 in data_dict[subjectID]
    }
    subjectIDs_timepoint2 = {
        subjectID for subjectID in data_dict if timepoint_2 in data_dict[subjectID]
    }

    # Find subjectIDs present only in one timepoint
    subjectIDs_only_in_timepoint1 = subjectIDs_timepoint1 - subjectIDs_timepoint2
    subjectIDs_only_in_timepoint2 = subjectIDs_timepoint2 - subjectIDs_timepoint1

    # Log warnings for subjectIDs present only in one timepoint
    if subjectIDs_only_in_timepoint1:
        logging.warning(
            f"The following subjectIDs are present only in timepoint '{timepoint_1}' and not in timepoint '{timepoint_2}': {subjectIDs_only_in_timepoint1}"
        )

    if subjectIDs_only_in_timepoint2:
        logging.warning(
            f"The following subjectIDs are present only in timepoint '{timepoint_2}' and not in timepoint '{timepoint_1}': {subjectIDs_only_in_timepoint2}"
        )

    results = []

    # Iterate over subjectIDs present in both timepoints
    common_subjectIDs = [
        subjectID
        for subjectID in data_dict
        if timepoint_1 in data_dict[subjectID] and timepoint_2 in data_dict[subjectID]
    ]

    if not common_subjectIDs:
        raise ValueError(
            f"No common subjectIDs found between the {timepoint_1} and {timepoint_2}."
        )
        sys.exit(42)

    for subjectID in common_subjectIDs:
        # Get DataFrames for each timepoint
        df_timepoint1 = data_dict[subjectID][timepoint_1]
        df_timepoint2 = data_dict[subjectID][timepoint_2]

        # Merge on contig and position
        # In pandas when both key columns ('gene_id' here) contain rows where the key is a null value, those rows will be matched against each other
        """
        df1 = pd.DataFrame({'key1': [1, 2, None], 'value1': ['A', 'B', 'C']})
        df2 = pd.DataFrame({'key1': [1, None, 3], 'value1': ['X', 'Y', 'Z']})

        merged_df = pd.merge(df1, df2, on='key1', how='inner')
        print(merged_df)
           key1 value1_x value1_y
        0   1.0        A        X
        1   NaN        C        Y
        """
        merged_df = pd.merge(
            df_timepoint1,
            df_timepoint2,
            on=["subjectID", "contig", "gene_id", "position", "replicate", "group"],
            suffixes=(f"_{timepoint_1}", f"_{timepoint_2}"),
            how="inner",
        )

        if merged_df.empty:
            logging.warning(f"No matching positions found for subjectID {subjectID}.")
            continue

        # Compute differences
        # Since the merge is "inner", only the positions common to both timepoints are present
        for nuc in NUCLEOTIDES:
            merged_df[f"{nuc}_diff"] = (
                merged_df[f"{nuc}_{timepoint_2}"] - merged_df[f"{nuc}_{timepoint_1}"]
            )

        # Select relevant columns
        columns_to_keep = (
            ["subjectID", "gene_id", "contig", "position", "replicate", "group"]
            + [f"{nuc}_{timepoint_1}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_{timepoint_2}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_diff" for nuc in NUCLEOTIDES]
        )

        results.append(merged_df[columns_to_keep])

    if not results:
        logging.error("No allele frequency changes calculated.")
        sys.exit(42)

    allele_changes = pd.concat(results, ignore_index=True)

    # Save the allele frequency changes
    # allele_changes.to_csv(
    #     os.path.join(output_dir, f"{mag_id}_allele_frequency_changes.tsv.gz"),
    #     sep="\t",
    #     index=False,
    #     compression="gzip",
    # )

    logging.info(
        f"Allele frequency changes saved to {output_dir}/{mag_id}_allele_frequency_changes.tsv.gz"
    )

    return allele_changes


def filter_zero_diff_positions(allele_changes, output_dir, mag_id):
    """
    Filters out positions where the sum of the difference values across all samples for all nucleotides is zero.

    Args:
        allele_changes (pd.DataFrame): DataFrame containing allele frequency changes with columns for contig, position,
                                       and nucleotide frequency differences (e.g., 'A_diff', 'T_diff', 'G_diff', 'C_diff').

    Returns:
        pd.DataFrame: A filtered DataFrame with positions where the sum of the difference values across all samples
                      for all nucleotides is zero removed.
    """
    logging.info(
        f"Identifying positions where sum of the difference values across all samples for all nucleotides is zero, called zero-diff positions"
    )

    diff_cols = [f"{nuc}_diff" for nuc in NUCLEOTIDES]

    # Group by (contig, position) and sum the _diff columns for each group.
    # This gives us, for each (contig, position), the total difference across all samples.
    grouped_sums = allele_changes.groupby(["contig", "position"], dropna=False)[
        diff_cols
    ].sum()

    # Find positions where ALL of the diff-sums are zero.
    is_all_zero = (
        (grouped_sums["A_frequency_diff"] == 0)
        & (grouped_sums["T_frequency_diff"] == 0)
        & (grouped_sums["G_frequency_diff"] == 0)
        & (grouped_sums["C_frequency_diff"] == 0)
    )
    # multi-index of (contig, gene_id, position)
    zero_positions = is_all_zero[is_all_zero].index

    logging.info(f"Found {len(zero_positions):,} zero-diff positions.")

    # Filter those positions OUT of allele_changes
    logging.info(f"Filtering zero-diff positions")
    ac_indexed = allele_changes.set_index(["contig", "position"], drop=False)
    keep_mask_ac = ~ac_indexed.index.isin(zero_positions)
    filtered_allele_changes = ac_indexed[keep_mask_ac].copy()
    filtered_allele_changes.reset_index(drop=True, inplace=True)

    # Save the allele frequency changes
    filtered_allele_changes.to_csv(
        os.path.join(
            output_dir, f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
        ),
        sep="\t",
        index=False,
        compression="gzip",
    )

    logging.info(
        f"Allele frequency changes with no zero diff positions saved to {output_dir}/{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
    )

    return filtered_allele_changes


def get_mean_change(allele_changes, mag_id, output_dir):
    """
    Calculate the mean changes in allele frequencies for subjectIDs present in the same replicate and group.

    Parameters:
        allele_changes (pd.DataFrame): DataFrame containing allele change information with columns
                                       ['contig', 'position', 'replicate', 'group', 'subjectID', 'gene_id' '<nuc>_diff'].
        mag_id (str): Identifier for the metagenome-assembled genome (MAG).
        output_dir (str): Directory where the output file will be saved.

    Returns:
        pd.DataFrame: DataFrame with mean changes in allele frequencies and count of unique subject IDs
                      for each group, with columns ['contig', 'position', 'replicate', 'group', 'gene_id',
                      '<nuc>_diff_mean', 'subjectID_count'].

    Notes:
        - The function groups the input DataFrame by ['contig', 'position', 'replicate', 'group', 'gene_id'] and
          calculates the mean of nucleotide differences and the count of unique subject IDs in each group.
        - The resulting DataFrame is saved as a compressed TSV file in the specified output directory.
    """
    logging.info(
        "Calculating mean changes in allele frequencies for subjectIDs present in the same replicate and group."
    )
    # Prepare the aggregation dictionary
    # https://pandas.pydata.org/docs/reference/api/pandas.core.groupby.DataFrameGroupBy.agg.html
    agg_dict = {f"{nuc}_diff": "mean" for nuc in NUCLEOTIDES}
    agg_dict["subjectID"] = "nunique"  # Count of unique subject IDs in each group

    # Perform the groupby and aggregation
    mean_changes_df = (
        allele_changes.groupby(
            ["contig", "gene_id", "position", "replicate", "group"], dropna=False
        )
        .agg(agg_dict)
        .reset_index()
    )

    # Rename the nucleotide columns to reflect mean changes
    mean_changes_df.rename(
        columns={f"{nuc}_diff": f"{nuc}_diff_mean" for nuc in NUCLEOTIDES}, inplace=True
    )

    # Optionally, rename 'subjectID' to 'subjectID_count' to make it clear it's a count
    mean_changes_df.rename(columns={"subjectID": "subjectID_count"}, inplace=True)

    mean_changes_df.to_csv(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_changes_mean.tsv.gz"),
        index=False,
        sep="\t",
        compression="gzip",
    )
    logging.info(
        f"Mean change in allele frequency for MAG {mag_id} saved to {output_dir}"
    )
    return mean_changes_df


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Analyze allele frequency and perform significance tests.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--magID",
        help="MAG ID to process",
        type=str,
        required=True,
        metavar="str",
    )
    parser.add_argument(
        "--mag_metadata_file",
        help="Path to metadata file",
        type=str,
        required=True,
        metavar="filepath",
    )
    parser.add_argument(
        "--fasta",
        help="Path to FASTA file with contigs",
        type=str,
        required=True,
        metavar="filepath",
    )
    parser.add_argument(
        "--breath_threshold",
        help="Breath threshold to use for MAGs.",
        type=float,
        default=0.1,
        metavar="float",
    )

    parser.add_argument(
        "--disable_zero_diff_filtering",
        help="Do not remove positions where change in allele frequency for each nucleotide's for in samples sums to zero across. Default is to do the filtering.",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--cpus",
        help=f"Number of processors to use.",
        default=cpu_count(),
        metavar="int",
        type=int,
    )

    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
        type=str,
        required=True,
        metavar="Directory path",
    )

    args = parser.parse_args()
    start_time = time.time()
    # Calculate MAG sizes
    mag_size_dict = calculate_mag_sizes(args.fasta)
    mag_id = args.magID

    # Load per-MAG metadata file and get required data structures
    metadata_dict, sample_files_with_mag_id = load_mag_metadata_file(
        args.mag_metadata_file, mag_id, args.breath_threshold
    )

    # Process the samples
    number_of_processes = min(args.cpus, len(sample_files_with_mag_id))

    logging.info(
        f"Processing {len(sample_files_with_mag_id)} samples for MAG {mag_id} using {number_of_processes} processes."
    )

    # Load sample data in parallel
    with Pool(
        processes=number_of_processes,
        initializer=init_worker,
        initargs=(metadata_dict, mag_size_dict),
    ) as pool:
        data_list = list(
            pool.imap_unordered(process_mag_files, sample_files_with_mag_id)
        )
    # Filter out None values (samples with breadth < 50%)
    data_list = [df for df in data_list if df is not None]

    if not data_list:
        logging.error(
            f"No samples for MAG {mag_id} passed the breadth threshold. Exiting...."
        )
        sys.exit(0)  # Exit the program

    data_dict = create_data_dict(data_list)

    # Release memory from data_list
    del data_list
    gc.collect()

    # save_nucleotide_frequencies(data_dict, args.output_dir, mag_id)

    allele_changes = calculate_allele_frequency_changes(
        data_dict, args.output_dir, mag_id
    )

    if not args.disable_zero_diff_filtering:
        logging.info("Filtering zero-diff positions.")
        allele_changes = filter_zero_diff_positions(
            allele_changes, args.output_dir, mag_id
        )
    else:
        logging.info("User disabled zero-diff position filtering.")

    get_mean_change(allele_changes, mag_id, args.output_dir)

    end_time = time.time()
    logging.info(f"Total time taken: {end_time-start_time:.2f} seconds")


if __name__ == "__main__":
    main()
