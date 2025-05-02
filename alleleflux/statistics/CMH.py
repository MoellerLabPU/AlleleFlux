#!/usr/bin/env python
"""
CMH Test Runner for AlleleFlux
-----------------------------
Performs Cochran-Mantel-Haenszel (CMH) tests on allele count data for metagenome-assembled genomes (MAGs),
leveraging R's mantelhaen.test via rpy2. Supports multiprocessing and progress tracking.

Usage:
    python cmh_test.py --input_df <input.tsv> --preprocessed_df <preprocessed.tsv> --mag_id <MAG> --output_dir <dir>
"""
import argparse
import logging
import os
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
from tqdm import tqdm

NUCLEOTIDES: List[str] = ["A", "T", "G", "C"]

# Rpy2 Setup. Enable pandas <-> R DataFrame conversion
pandas2ri.activate()


def load_and_filter_data(
    input_df_path: str, preprocessed_df_path: str, mag_id: str, dtype_map: dict
) -> pd.DataFrame:
    """
    Load raw allele count data and filter it to include only positions present in a preprocessed dataset.

    Parameters:
        input_df_path (str): Path to the raw allele count data file (tab-separated values).
        preprocessed_df_path (str): Path to the preprocessed positions file (tab-separated values).
        mag_id (str): Identifier for the metagenome-assembled genome (MAG), used for logging and error messages.

    Returns:
        pd.DataFrame: Filtered DataFrame containing only rows from the raw count data whose (contig, position)
                      pairs are present in the preprocessed positions.

    Raises:
        ValueError: If no positions remain after filtering by the preprocessed data.
    """
    # Read the raw allele count data
    logging.info(f"Loading raw count data from {input_df_path}")
    raw_counts_df = pd.read_csv(
        input_df_path,
        sep="\t",
        usecols=dtype_map.keys(),
        dtype=dtype_map,
        index_col=["contig", "position"],
        memory_map=True,
        # low_memory=False,
        # nrows=20000000,
    )
    # Log basic row count quickly
    logging.info(f"Loaded {raw_counts_df.shape[0]:,} rows of raw count data.")
    # Detailed stats are debug level for performance
    logging.debug(
        f"Distinct contigs: {raw_counts_df.index.get_level_values('contig').nunique():,}",
    )
    logging.debug(
        f"Distinct positions: {raw_counts_df.index.nunique():,}",
    )

    # Read the preprocessed positions for filtering
    logging.info(f"Loading preprocessed positions from {preprocessed_df_path}")
    preprocessed_df = pd.read_csv(
        preprocessed_df_path,
        sep="\t",
        usecols=["contig", "position"],
    )
    logging.info(f"Loaded {preprocessed_df.shape[0]:,} preprocessed positions.")

    # Get unique (contig, position) pairs to keep
    valid_positions = preprocessed_df.drop_duplicates(subset=["contig", "position"])
    valid_positions_index = pd.MultiIndex.from_frame(valid_positions)

    # Filter Raw Counts by Preprocessed Positions
    logging.info("Filtering raw counts to include only preprocessed positions.")

    # Keep only rows whose (contig, position) index is in the valid_positions_index
    filtered_counts_df = raw_counts_df[
        raw_counts_df.index.isin(valid_positions_index)
    ].reset_index()

    logging.info(
        f"Filtered to {filtered_counts_df.shape[0]:,} rows after applying preprocessed positions."
    )

    if filtered_counts_df.empty:
        raise ValueError(
            f"No positions remaining after filtering by preprocessed data for MAG {mag_id}. Exiting."
        )

    logging.debug(
        f"{len(filtered_counts_df['contig'].unique()):,} contigs and "
        f"{len(filtered_counts_df.groupby(['contig', 'position'])):,} unique positions remaining after filtering."
    )
    return filtered_counts_df


def process_group(
    args: Tuple[Tuple[str, str, int], pd.DataFrame],
    min_sample_num: int,
    group_1_name: str,
    group_2_name: str,
    mag_id: str,
    timepoint: Optional[str] = None,
) -> dict:
    """
    Processes a group of allele count data for a specific gene position and performs the Cochran-Mantel-Haenszel (CMH) test.

    Parameters:
        args (Tuple[Tuple[str, str, int], pd.DataFrame]):
            A tuple containing:
                - name (Tuple[str, str, int]): Identifiers for the group (contig, gene_id, position).
                - df_group (pd.DataFrame): DataFrame containing allele counts for the group.
        min_sample_num (int): Minimum number of samples required for the CMH test.
        group_1_name (str): Name of the first comparison group.
        group_2_name (str): Name of the second comparison group.
        mag_id (str): Identifier for the metagenome-assembled genome (MAG).
        timepoint (Optional[str]): The timepoint identifier for longitudinal data.

    Returns:
        dict or None:
            - A dictionary with test results and metadata if the group passes filtering and the CMH test runs successfully.
            - Returns None if there is insufficient data for the test.
            - If an error occurs during the CMH test, returns a dictionary with NaN p-value and error notes.
    """
    name, df_group = args  # Unpack the args tuple
    contig, gene_id, position = name  # Extract identifiers

    # Prepare the input for CMH test
    pivoted_filtered = prepare_cmh_input(
        df_group, group_1_name, group_2_name, min_sample_num
    )
    if pivoted_filtered is None:
        return None  # Skip if not enough data

    # logging.info(f"Running CMH for {contig}:{position}")
    # Run the CMH test in R
    try:
        cmh_result = run_cmh_test_in_r(pivoted_filtered)
    except RRuntimeError as e:
        logging.error(
            f"Error running CMH test for {contig}:{position} in MAG {mag_id}: {e}"
        )
        # Return NaN p-value with error note
        result = {
            "mag_id": mag_id,
            "contig": contig,
            "gene_id": gene_id,
            "position": position,
            "num_pairs": len(pivoted_filtered),
            "p_value_CMH": np.nan,
            "notes": f"RRuntimeError: {str(e).rstrip('\n')}",
        }
        if timepoint is not None:
            result["time"] = timepoint
        return result

    result = {
        "mag_id": mag_id,
        "contig": contig,
        "gene_id": gene_id,
        "position": position,
        "num_pairs": len(pivoted_filtered),
        "p_value_CMH": cmh_result,
        "notes": "",
    }
    if timepoint is not None:
        result["time"] = timepoint
    return result


def get_informative_strata(
    pivoted_df: pd.DataFrame, group_1_name: str, group_2_name: str
) -> List[str]:
    """
    Identifies nucleotide bases that are informative (i.e., present in at least one of two groups) within a pivoted DataFrame.

    Args:
        pivoted_df (pd.DataFrame): A DataFrame where columns are named as '{group}_{base}' for each group and nucleotide base.
        group_1_name (str): The name of the first group.
        group_2_name (str): The name of the second group.

    Returns:
        List[str]: A list of nucleotide bases that are present (sum > 0) in at least one of the two groups.
    """
    informative: List[str] = []
    for base in NUCLEOTIDES:
        # Build column names for each group
        col1 = f"{group_1_name}_{base}"
        col2 = f"{group_2_name}_{base}"
        # Only keep bases present in both groups
        if col1 in pivoted_df.columns and col2 in pivoted_df.columns:
            if pivoted_df[col1].sum() > 0 or pivoted_df[col2].sum() > 0:
                informative.append(base)
    return informative


def prepare_cmh_input(
    df_group: pd.DataFrame,
    group_1_name: str,
    group_2_name: str,
    min_sample_num: int,
) -> Optional[pd.DataFrame]:
    """
    Prepares input data for the Cochran-Mantel-Haenszel (CMH) test by aggregating nucleotide counts
    for two specified groups and filtering for informative bases and shared replicates.

    The function performs the following steps:
    - Splits the input DataFrame into two groups based on the provided group names.
    - Aggregates nucleotide counts by contig, position, group, and replicate.
    - Identifies replicates shared between both groups and filters out those not present in both.
    - Pivots the data so that each row corresponds to a replicate, and columns represent nucleotide counts per group.
    - Identifies informative nucleotide bases (those with nonzero counts in at least one group).
    - Filters out replicates with total nucleotide count less than or equal to 1.
    - Ensures that the number of replicates meets the minimum sample requirement.

    Parameters
    ----------
    df_group : pd.DataFrame
        Input DataFrame containing columns: 'contig', 'position', 'group', 'replicate', and nucleotide counts.
    group_1_name : str
        Name of the first group to compare.
    group_2_name : str
        Name of the second group to compare.
    min_sample_num : int
        Minimum number of shared replicates required to proceed.

    Returns
    -------
    Optional[pd.DataFrame]
        A DataFrame formatted for CMH test input, containing only informative bases and shared replicates,
        or None if requirements are not met (e.g., insufficient replicates or informative bases).
    """
    # Split group into two based on group name
    group1_df = df_group[df_group["group"] == group_1_name]
    group2_df = df_group[df_group["group"] == group_2_name]

    if group1_df.empty or group2_df.empty:
        return None  # If either group is missing, cannot proceed

    # Combine both groups for aggregation
    combined_df = pd.concat([group1_df, group2_df])

    # Aggregate nucleotide counts by contig, position, group, and replicate
    # Get the sum of each nucleotide for each replicate
    aggregated_df = (
        combined_df.groupby(["contig", "position", "group", "replicate"], dropna=False)[
            NUCLEOTIDES
        ]
        .sum()
        .reset_index()
    )

    # Find unique replicates present in both groups
    reps_by_group = aggregated_df.groupby("group")["replicate"].unique()

    # Get shared replicates
    shared_reps = np.intersect1d(
        reps_by_group[group_1_name], reps_by_group[group_2_name]
    )

    # Skip if not enough replicates
    if len(shared_reps) < min_sample_num:
        return None

    # Keep only replicates present in both groups
    subset = aggregated_df[aggregated_df["replicate"].isin(shared_reps)]

    # Pivot so each row is a replicate, columns are group_nucleotide (e.g., control_A)
    pivoted = (
        subset.pivot(
            index=["contig", "position", "replicate"],
            columns="group",
            values=NUCLEOTIDES,
        )
        .fillna(0)
        .astype(int)
    )
    # Flatten MultiIndex columns: (nucleotide, group) -> 'group_nucleotide'
    pivoted.columns = [f"{group}_{nt}" for nt, group in pivoted.columns]
    pivoted = pivoted.reset_index()

    # Keeps the base if at least one of the groups have a nonzero sum.
    informative_bases = get_informative_strata(pivoted, group_1_name, group_2_name)
    # P-value is NaN when only 2 bases are present and one of them is 0
    if len(informative_bases) < 2:
        return None  # Need at least two informative bases for CMH

    # Build list of columns to keep: contig, position, replicate, and informative group_nucleotide columns
    keep_cols = ["contig", "position", "replicate"]
    for base in informative_bases:
        keep_cols += [f"{group_1_name}_{base}", f"{group_2_name}_{base}"]
    pivoted_subset = pivoted[keep_cols].copy()

    nucleotide_cols = [
        col
        for col in pivoted_subset.columns
        if any(col.endswith(f"_{nuc}") for nuc in NUCLEOTIDES)
    ]
    # Sum all nucleotide columns to get the total count for each replicate
    pivoted_subset["total"] = pivoted_subset[nucleotide_cols].sum(axis=1)

    # Filter out replicates with total <= 1, else we get the error:
    # sample size in each stratum must be > 1
    pivoted_filtered = (
        pivoted_subset[pivoted_subset["total"] > 1].copy().drop(columns=["total"])
    )

    # Skip if not enough replicates after filtering
    if len(pivoted_filtered) < min_sample_num:
        return None

    return pivoted_filtered


def run_cmh_test_in_r(pivoted_filtered: pd.DataFrame) -> Optional[float]:
    """
    Runs the Cochran-Mantel-Haenszel (CMH) test in R on a given pivoted and filtered pandas DataFrame.

    This function transfers the provided DataFrame to the R environment, reshapes it for analysis,
    constructs a 3D contingency table, and performs the CMH test using R's `mantelhaen.test`.
    The DataFrame is expected to have columns: 'contig', 'position', 'replicate', and additional columns
    with group and nucleotide information separated by an underscore (e.g., 'fat_A', 'control_A').

    Args:
        pivoted_filtered (pd.DataFrame): A pandas DataFrame containing the data to be analyzed,
            with specific columns as described above.

    Returns:
        Optional[float]: The p-value from the CMH test if successful, otherwise None.
    """
    # Transfer DataFrame to R
    ro.globalenv["r_df"] = pivoted_filtered
    # R code is generic: group names are parsed from column names
    r(
        f"""
        library(tidyr)
        # Convert wide to long format for R
        melted <- pivot_longer(
            r_df,
            cols = -c(contig, position, replicate),
            names_to = c("group", "nucleotide"),
            names_sep = "_",
            values_to = "count"
        )
        melted$count <- as.integer(melted$count)
        # Build 3D contingency table
        # Last dimension needs to be the 'strata' ie. replicate. Firt 2 columns are interchangeable
        tab3d <- xtabs(count ~ group + nucleotide + replicate, data = melted)
        # Run CMH test
        cmh_result <- mantelhaen.test(tab3d, alternative = c("two.sided"), correct=FALSE)$p.value
        """
    )
    cmh_result = ro.globalenv["cmh_result"]
    return cmh_result[0]


def run_cmh_tests(
    df: pd.DataFrame,
    mag_id: str,
    cpus: int,
    min_sample_num: int,
    timepoint: Optional[str] = None,
) -> Optional[pd.DataFrame]:
    """
    Run CMH test workflow on the given DataFrame.

    Args:
        df: DataFrame containing the data to process
        mag_id: MAG identifier
        cpus: Number of CPUs to use
        min_sample_num: Minimum number of samples required
        timepoint: Timepoint identifier for longitudinal data

    Returns:
        Optional[pd.DataFrame]: DataFrame with results if successful, None otherwise
    """
    groups = df["group"].unique()
    if len(groups) != 2:
        raise ValueError(
            f"Expected exactly 2 groups, but found {len(groups)} groups: {groups}"
        )
    group_1, group_2 = groups
    logging.info(f"Grouping data by contig, gene_id and position")
    grouped = df.groupby(["contig", "gene_id", "position"], dropna=False)
    count = len(grouped)
    logging.info(
        f"Analyzing {count:,} positions for CMH tests between {group_1} and {group_2}{f' (timepoint: {timepoint})' if timepoint else ''} using {cpus} cores"
    )
    with Pool(processes=cpus) as pool:
        worker = partial(
            process_group,
            min_sample_num=min_sample_num,
            group_1_name=group_1,
            group_2_name=group_2,
            mag_id=mag_id,
            timepoint=timepoint,
        )
        results = [
            res
            for res in tqdm(
                pool.imap_unordered(worker, grouped),
                total=count,
                desc=f"Running CMH test{(f' {timepoint}' if timepoint else '')}",
            )
            if res is not None
        ]
    if not results:
        logging.warning(
            f"No results for {mag_id} {f' (timepoint: {timepoint})' if timepoint else ''}."
        )
        return None

    df_out = pd.DataFrame(results)
    return df_out


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s %(filename)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Run Cochran-Mantel-Haenszel test for a MAG across different samples using R's mantelhaen.test via rpy2. Uses raw counts and filters by preprocessed positions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_df",
        required=True,
        help="Path to input allele frequency dataframe. For single timepoint data, preprocessed_df from preprocess_two_sample.py can be used for --input_df.",
        type=str,
    )
    parser.add_argument(
        "--preprocessed_df",
        help="Path to the filtered dataframe from preprocess_two_sample.py, used to filter positions. Only required for longitudinal data. For single timepoint data, preprocessed_df can be used for --input_df.",
        type=str,
    )
    parser.add_argument(
        "--min_sample_num",
        type=int,
        default=4,
        help="Minimum number of strata (replicates) required",
    )
    parser.add_argument(
        "--mag_id",
        help="MAG ID to process",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--data_type",
        help="Type of data to analyze: longitudinal or single",
        type=str,
        choices=["longitudinal", "single"],
        default="longitudinal",
    )
    parser.add_argument(
        "--cpus",
        help="Number of processors to use.",
        default=cpu_count(),
        type=int,
    )
    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
        type=str,
        required=True,
    )
    args = parser.parse_args()

    # Define data types for columns
    dtype_map = {
        "contig": str,
        "gene_id": str,
        "position": int,
        "group": str,
        "replicate": str,
        **{nuc: int for nuc in NUCLEOTIDES},
    }

    # Load and possibly filter data
    if args.data_type == "single":
        logging.info(
            "Datatype set to single. Loading input dataframe without any filtering."
        )
        df = pd.read_csv(
            args.input_df,
            sep="\t",
            usecols=dtype_map.keys(),
            dtype=dtype_map,
        )

        run_cmh_tests(
            df=df,
            mag_id=args.mag_id,
            cpus=args.cpus,
            min_sample_num=args.min_sample_num,
        )
    else:
        # Longitudinal
        dtype_map["time"] = str
        if args.preprocessed_df:
            logging.info(
                "Datatype set to longitudinal. Loading input dataframe with filtering."
            )
            df = load_and_filter_data(
                args.input_df, args.preprocessed_df, args.mag_id, dtype_map
            )
        else:
            logging.info(
                "Datatype set to longitudinal. No preprocessed dataframe provided. Using input_df directly."
            )
            df = pd.read_csv(
                args.input_df,
                sep="\t",
                usecols=dtype_map.keys(),
                dtype=dtype_map,
            )

        timepoints = df["time"].unique()
        if len(timepoints) != 2:
            raise ValueError(
                f"Expected exactly 2 unique timepoints, but found {len(timepoints)}: {timepoints}"
            )

        # Process each timepoint and collect results
        all_results = []
        for tp in timepoints:
            logging.info(f"Processing timepoint: {tp}")
            df_tp = df[df["time"] == tp].copy()

            # For individual timepoint files
            tp_results = run_cmh_tests(
                df=df_tp,
                mag_id=args.mag_id,
                cpus=args.cpus,
                min_sample_num=args.min_sample_num,
                timepoint=tp,
            )

            if tp_results is not None:
                all_results.append(tp_results)

        # Combine results from all timepoints into a single dataframe
        if all_results:
            os.makedirs(args.output_dir, exist_ok=True)
            combined_results = pd.concat(all_results, ignore_index=True)
            combined_out_file = os.path.join(
                args.output_dir, f"{args.mag_id}_cmh.tsv.gz"
            )
            logging.info(f"Saving combined CMH results to {combined_out_file}")
            combined_results.to_csv(
                combined_out_file, sep="\t", index=False, compression="gzip"
            )
            logging.info(f"Saved combined CMH test results to {combined_out_file}")
        else:
            logging.warning(f"No results from any timepoint for {args.mag_id}.")


if __name__ == "__main__":
    main()
