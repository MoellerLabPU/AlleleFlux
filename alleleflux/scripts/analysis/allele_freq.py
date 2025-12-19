import argparse
import gc
import logging
import os
import sys
import time
from multiprocessing import Pool, cpu_count

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]

logger = logging.getLogger(__name__)


def load_qc_results(qc_file_path, mag_id):
    """
    Load QC results and return samples that passed coverage threshold.

    Parameters:
        qc_file_path (str): Path to the QC TSV file for this MAG
        mag_id (str): MAG identifier (for logging)

    Returns:
        pd.DataFrame: Filtered DataFrame containing only samples that passed coverage threshold,
                     with columns: sample_id, file_path, group, subjectID, replicate, time (if available),
                     breadth, average_coverage, genome_size, and other QC metrics

    Raises:
        ValueError: If no samples passed the coverage threshold
        FileNotFoundError: If QC file doesn't exist
    """
    logger.info(f"Loading QC results from {qc_file_path}")

    qc_df = pd.read_csv(qc_file_path, sep="\t", dtype={"sample_id": str})

    # Filter to samples that passed coverage threshold
    if "coverage_threshold_passed" not in qc_df.columns:
        raise ValueError(
            f"QC file {qc_file_path} missing 'coverage_threshold_passed' column. "
            "Ensure quality_control.py was run properly."
        )

    passed_df = qc_df[qc_df["coverage_threshold_passed"] == True].copy()

    if passed_df.empty:
        raise ValueError(
            f"No samples passed coverage threshold for MAG {mag_id}. "
            f"Total samples in QC file: {len(qc_df)}"
        )

    logger.info(
        f"Loaded {len(passed_df)} samples (out of {len(qc_df)} total) that passed coverage and breadth thresholds"
    )

    return passed_df


def init_worker(metadata, data_type):
    """
    Initialize worker process with metadata.

    This function sets up global dictionary for metadata
    that can be accessed by worker processes.

    Parameters:
        metadata (dict): A dictionary containing metadata information with QC metrics.
        data_type (str): Type of data analysis ("single" or "longitudinal").

    Returns:
        None
    """
    global metadata_dict
    global DATA_TYPE
    metadata_dict = metadata
    DATA_TYPE = data_type


def process_mag_files(args):
    """
    Processes a MAG profile file using metadata from QC results.

    This function reads a coverage profile and adds metadata from QC analysis.
    Since QC has already verified breadth and coverage thresholds, this function
    simply loads the profile and enriches it with metadata.

    Parameters:
        args (tuple): A tuple containing:
            - sample_id (str): The sample identifier
            - filepath (str): Path to the coverage profile file
            - mag_id (str): MAG identifier

    Returns:
        pd.DataFrame or None: DataFrame with nucleotide frequencies and metadata if successful,
                             None if file cannot be read

    Notes:
        - Reads the profile file from filepath
        - Adds metadata columns from the global metadata_dict (from QC results)
        - Inserts MAG_ID as the first column
        - Calculates nucleotide frequencies
        - All samples have already passed breadth and coverage thresholds (verified in QC)
    """
    sample_id, filepath, mag_id = args

    df = pd.read_csv(filepath, sep="\t", dtype={"gene_id": str})

    # Add sample_id column
    df["sample_id"] = sample_id

    # Add metadata columns from QC results
    metadata_info = metadata_dict.get(sample_id)
    if metadata_info is None:
        logger.error(f"Sample ID not found for sample {sample_id}")
        return None

    df["group"] = metadata_info["group"]
    df["subjectID"] = metadata_info["subjectID"]
    df["replicate"] = metadata_info["replicate"]

    # Add QC metrics from metadata
    df["breadth"] = metadata_info["breadth"]
    df["genome_size"] = metadata_info["genome_size"]
    df["average_coverage"] = metadata_info["average_coverage"]

    # For longitudinal data, add time column if available
    if DATA_TYPE == "longitudinal" and "time" in metadata_info:
        df["time"] = metadata_info["time"]

    # Insert MAG_ID as the first column
    df.insert(0, "MAG_ID", mag_id)

    # Calculate nucleotide frequencies
    df = calculate_frequencies(df)
    return df


def calculate_frequencies(df):
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

    # Calculate frequencies directly
    df["A_frequency"] = df["A"] / total_coverage
    df["T_frequency"] = df["T"] / total_coverage
    df["G_frequency"] = df["G"] / total_coverage
    df["C_frequency"] = df["C"] / total_coverage

    return df


def save_allele_frequencies(data_dict, output_dir, mag_id):
    mag_df = pd.concat(
        [df for subject_dict in data_dict.values() for df in subject_dict.values()],
        ignore_index=True,
    )

    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Saving nucleotide frequencies for MAG {mag_id} to {output_dir}")
    mag_df.to_csv(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_longitudinal.tsv.gz"),
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
    logger.info("Identifying unique timepoints.")

    unique_timepoints = set()
    for subject_data in data_dict.values():
        unique_timepoints.update(subject_data.keys())
    if len(unique_timepoints) != 2:
        raise ValueError(
            f"Expected exactly 2 unique timepoints, found {len(unique_timepoints)}."
        )

    # Unpack the two timepoints
    timepoint_1, timepoint_2 = unique_timepoints

    logger.info(
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
        logger.warning(
            f"The following subjectIDs are present only in timepoint '{timepoint_1}' and not in timepoint '{timepoint_2}': {subjectIDs_only_in_timepoint1}"
        )

    if subjectIDs_only_in_timepoint2:
        logger.warning(
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
    logger.info(f"Common subjectIDs: {common_subjectIDs}")

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
            logger.warning(f"No matching positions found for subjectID {subjectID}.")
            continue

        # Compute differences
        # Since the merge is "inner", only the positions common to both timepoints are present
        for nuc in NUCLEOTIDES:
            merged_df[f"{nuc}_diff"] = (
                merged_df[f"{nuc}_{timepoint_2}"] - merged_df[f"{nuc}_{timepoint_1}"]
            )

        # Calculate combined total coverage
        merged_df["total_coverage_combined"] = (
            merged_df[f"total_coverage_{timepoint_1}"]
            + merged_df[f"total_coverage_{timepoint_2}"]
        )

        # Select relevant columns
        columns_to_keep = (
            ["subjectID", "gene_id", "contig", "position", "replicate", "group"]
            + [
                f"total_coverage_{timepoint_1}",
                f"total_coverage_{timepoint_2}",
                "total_coverage_combined",
            ]
            + [f"{nuc}_{timepoint_1}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_{timepoint_2}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_diff" for nuc in NUCLEOTIDES]
        )

        results.append(merged_df[columns_to_keep])

    if not results:
        logger.error("No allele frequency changes calculated.")
        sys.exit(42)

    allele_changes = pd.concat(results, ignore_index=True)

    # Save the allele frequency changes
    allele_changes.to_csv(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_changes.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
    )

    logger.info(
        f"Allele frequency changes saved to {output_dir}/{mag_id}_allele_frequency_changes.tsv.gz"
    )

    return allele_changes


def filter_constant_positions(allele_df, output_dir, mag_id, data_type):
    """
    Filter positions where allele frequencies do not change across samples.

    For single data type, it removes positions where all nucleotide frequencies
    are constant across all samples.

    For longitudinal data type, it removes 'zero-diff positions' where the sum
    of frequency differences across all samples for all nucleotides is zero.

    Parameters
    ----------
    allele_df : pandas.DataFrame
        DataFrame containing allele frequency data with columns for nucleotides (A, T, G, C)
        and positions information (contig, position)
    output_dir : str
        Directory where filtered output files will be saved
    mag_id : str
        Identifier for the metagenome-assembled genome (MAG)
    data_type : str
        Type of analysis to perform, either 'single' or 'longitudinal'

    Returns
    -------
    pandas.DataFrame or None
        For 'longitudinal' data type, returns DataFrame with zero-diff positions filtered out.
        For 'single' data type, returns None (results are saved to file only).

    Notes
    -----
    - For 'single' data type, saves filtered data to '[output_dir]/[mag_id]_allele_frequency_no_constant.tsv.gz'
    - For 'longitudinal' data type, saves filtered data to '[output_dir]/[mag_id]_allele_frequency_changes_no_zero-diff.tsv.gz'
    """

    grouped_df = allele_df.groupby(["contig", "position"], dropna=False)

    if data_type == "single":
        logger.info(
            f"Identifying positions where allele frequency values across all samples for each nucleotide is constant"
        )
        # Compute the number of unique values for each frequency column
        agg = grouped_df[NUCLEOTIDES].nunique()

        # Create a boolean mask that is True if at least one frequency column has more than one unique value
        groups_to_keep = (agg > 1).any(axis=1)

        positions_kept = groups_to_keep.sum()
        total_positions = agg.shape[0]
        positions_removed = total_positions - positions_kept

        logger.info(f"Found {positions_removed:,} constant positions. Filtering...")

        # groups_to_keep is a Series with MultiIndex (contig, position) and boolean values.
        # Select the groups (the index) where the condition is True.
        keep_index = groups_to_keep[groups_to_keep].index

        # Filter the original DataFrame to keep only rows that belong to the selected groups.
        filtered_df = (
            allele_df.set_index(["contig", "position"]).loc[keep_index].reset_index()
        )
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, Positions removed: {positions_removed:,}"
        )

        # Save the allele frequency changes
        filtered_df.to_csv(
            os.path.join(output_dir, f"{mag_id}_allele_frequency_no_constant.tsv.gz"),
            sep="\t",
            index=False,
            compression="gzip",
        )
        logger.info(
            f"Allele frequency changes with no constant positions saved to {output_dir}/{mag_id}_allele_frequency_no_constant.tsv.gz"
        )
        return None

    elif data_type == "longitudinal":

        logger.info(
            f"Identifying positions where sum of the difference values across all samples for all nucleotides is zero, called zero-diff positions"
        )

        diff_cols = [f"{nuc}_diff" for nuc in NUCLEOTIDES]

        # Group by (contig, position) and sum the _diff columns for each group.
        # This gives us, for each (contig, position), the total difference across all samples.
        grouped_sums = grouped_df[diff_cols].sum()

        # Find positions where ALL of the diff-sums are zero.
        is_all_zero = (
            (grouped_sums["A_frequency_diff"] == 0)
            & (grouped_sums["T_frequency_diff"] == 0)
            & (grouped_sums["G_frequency_diff"] == 0)
            & (grouped_sums["C_frequency_diff"] == 0)
        )
        # multi-index of (contig, gene_id, position)
        zero_positions = is_all_zero[is_all_zero].index

        logger.info(f"Found {len(zero_positions):,} zero-diff positions.")

        # Filter those positions OUT of allele_df
        logger.info(f"Filtering zero-diff positions")
        ac_indexed = allele_df.set_index(["contig", "position"], drop=False)
        keep_mask_ac = ~ac_indexed.index.isin(zero_positions)
        filtered_allele_changes = ac_indexed[keep_mask_ac].copy()
        filtered_allele_changes.reset_index(drop=True, inplace=True)

        total_positions = grouped_sums.shape[0]
        positions_removed = len(zero_positions)
        positions_kept = total_positions - positions_removed
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, Positions removed: {positions_removed:,}"
        )

        # Save the allele frequency changes
        filtered_allele_changes.to_csv(
            os.path.join(
                output_dir, f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
            ),
            sep="\t",
            index=False,
            compression="gzip",
        )

        logger.info(
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
    logger.info(
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
    logger.info(
        f"Mean change in allele frequency for MAG {mag_id} saved to {output_dir}"
    )
    return mean_changes_df


def process_single_data(data_list, output_dir, mag_id, disable_filtering):
    allele_df = pd.concat(data_list, ignore_index=True)

    output_fpath = os.path.join(output_dir, f"{mag_id}_allele_frequency_single.tsv.gz")
    logger.info(f"Writing allele frequencies (single data) to {output_fpath}")
    allele_df.to_csv(output_fpath, sep="\t", index=False, compression="gzip")

    if not disable_filtering:
        logger.info("Filtering constant allele frequency positions (single data).")

        # Call your vectorized filtering function for single data.
        filtered_df = filter_constant_positions(
            allele_df, output_dir, mag_id, data_type="single"
        )
    else:
        logger.info("User disabled filtering of constant positions.")


def process_longitudinal_data(data_list, output_dir, mag_id, disable_filtering):
    # Create a dictionary of DataFrames for each subjectID and timepoint
    data_dict = create_data_dict(data_list)

    # Free memory from data_list.
    del data_list
    gc.collect()

    output_fpath = os.path.join(
        output_dir, f"{mag_id}_allele_frequency_longitudinal.tsv.gz"
    )
    logger.info(f"Writing allele frequencies (longitudinal data) to {output_fpath}")
    save_allele_frequencies(data_dict, output_dir, mag_id)

    allele_changes = calculate_allele_frequency_changes(data_dict, output_dir, mag_id)

    if not disable_filtering:
        logger.info("Filtering zero-diff positions for longitudinal data.")
        allele_changes = filter_constant_positions(
            allele_changes, output_dir, mag_id, data_type="longitudinal"
        )
    else:
        logger.info("User disabled zero-diff filtering for longitudinal data.")

    get_mean_change(allele_changes, mag_id, output_dir)


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Analyze allele frequency using QC-filtered samples.",
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
        "--qc_file",
        help="Path to QC results TSV file (output from quality_control.py) for this MAG",
        type=str,
        required=True,
        metavar="filepath",
    )
    parser.add_argument(
        "--data_type",
        help="Is the data from a single timepoint or from a time series (longitudinal)?",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )

    parser.add_argument(
        "--disable_zero_diff_filtering",
        dest="disable_filtering",
        help="""
        For single: Do not remove positions where all nucleotide frequencies are constant across all samples.
        For longitudinal: Do not remove positions where change in allele frequency 
        for each nucleotide's samples sums to zero across. Default is to do the filtering.
        """,
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--cpus",
        help="Number of processors to use.",
        type=int,
        default=cpu_count(),
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

    mag_id = args.magID

    # Load QC results and filter to passing samples
    qc_df = load_qc_results(args.qc_file, mag_id)

    # Build metadata dict from QC results (sample_id -> metadata with QC metrics)
    metadata_dict = {}
    sample_file_tuples = []

    for _, row in qc_df.iterrows():
        sample_id = str(row["sample_id"])

        # Build metadata dict entry with all relevant info
        meta = {
            "group": row["group"],
            "subjectID": row["subjectID"],
            "replicate": row["replicate"],
            "breadth": row["breadth"],
            "genome_size": row["genome_size"],
            "average_coverage": row["average_coverage"],
        }

        # Add time if available (longitudinal data)
        if "time" in row and pd.notna(row["time"]):
            meta["time"] = row["time"]

        metadata_dict[sample_id] = meta

        # Build tuple for processing: (sample_id, file_path, mag_id)
        sample_file_tuples.append((sample_id, row["file_path"], mag_id))

    # Process the samples
    number_of_processes = min(args.cpus, len(sample_file_tuples))

    logger.info(
        f"Processing {len(sample_file_tuples)} samples for MAG {mag_id} using {number_of_processes} processes."
    )

    # Load sample data in parallel
    with Pool(
        processes=number_of_processes,
        initializer=init_worker,
        initargs=(metadata_dict, args.data_type),
    ) as pool:
        data_list = list(pool.imap_unordered(process_mag_files, sample_file_tuples))

    # Filter out None values (in case of read failures)
    data_list = [df for df in data_list if df is not None]

    if not data_list:
        logger.error(
            f"No sample profiles could be loaded for MAG {mag_id}. Check file paths in QC results."
        )
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    if args.data_type == "single":
        process_single_data(data_list, args.output_dir, mag_id, args.disable_filtering)
    elif args.data_type == "longitudinal":
        process_longitudinal_data(
            data_list, args.output_dir, mag_id, args.disable_filtering
        )

    end_time = time.time()
    logger.info(f"Total time taken: {end_time-start_time:.2f} seconds")


if __name__ == "__main__":
    main()
