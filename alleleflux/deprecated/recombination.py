#!/usr/bin/env python3

import argparse
import logging
import sys
import time
from functools import reduce

import pandas as pd
from utilities import extract_test_columns


def sliding_window_significance(
    df, test_columns_dict, window_size=100, p_value_threshold=0.05
):
    """
    Calculate the significance of tests in sliding windows across a genomic dataset.

    Parameters:
    df (pd.DataFrame): DataFrame containing genomic data with columns 'contig' and 'position'.
    test_columns_dict (dict): Dictionary where keys are test names and values are lists of column names containing p-values for those tests.
    window_size (int, optional): Size of the sliding window in base pairs. Default is 100.
    p_value_threshold (float, optional): Threshold for determining significance of p-values. Default is 0.05.

    Returns:
    pd.DataFrame: DataFrame with columns for each test's total sites, significant sites, and percentage of significant sites per window,
                  along with window start and end positions.
    """

    df = df.copy()  # do not mutate original

    # Create window_index from 0-based position
    df["window_index"] = df["position"] // window_size

    results = []
    for test_name, test_cols in test_columns_dict.items():
        subset_cols = ["contig", "position", "window_index"] + test_cols
        subdf = df[subset_cols].copy()

        # Drop rows that are all NaN in p-value columns (optional)
        subdf.dropna(subset=test_cols, how="all", inplace=True)

        # Check for NaNs in all p-value columns and drop
        if subdf[test_cols].isnull().values.any():
            raise ValueError("NaNs found in p-value column.")

        # Mark significant if ANY p-value < threshold
        significant_col = f"is_significant_{test_name}"
        subdf[significant_col] = subdf[test_cols].lt(p_value_threshold).any(axis=1)

        # Group by (contig, window_index)
        grouped = subdf.groupby(["contig", "window_index"], dropna=False)

        # Total sites = number of present sites per group
        total_sites = grouped.size()

        # Number of significant sites per group
        significant_sites = grouped[significant_col].sum()

        # Percentage of significant sites per group
        percentage = (significant_sites / total_sites) * 100

        test_df = pd.DataFrame(
            {
                "contig": total_sites.index.get_level_values("contig"),
                "window_index": total_sites.index.get_level_values("window_index"),
                f"total_sites_{test_name}": total_sites.values,
                f"significant_sites_{test_name}": significant_sites.values,
                f"score_{test_name} (%)": percentage.values,
            }
        )
        results.append(test_df)

    # Merge all test DataFrames
    final_df = reduce(
        lambda left, right: pd.merge(
            left, right, on=["contig", "window_index"], how="outer"
        ),
        results,
    )

    # Add window_start, window_end
    final_df["window_start"] = final_df["window_index"] * window_size
    final_df["window_end"] = final_df["window_start"] + window_size - 1

    return final_df


def filter_contigs(contigs_file, mag_id, length_cutoff):
    """
    Filters contigs from a given file based on MAG ID and length cutoff.

    Parameters:
        contigs_file (str): Path to the file containing contig lengths.
        mag_id (str): The MAG ID to filter contigs for.
        length_cutoff (int): The minimum length of contigs to retain.

    Returns:
        set: A set of valid contig names that meet the criteria.

    Raises:
        SystemExit: If the contigs_file does not have the required columns,
                    if no contigs are found for the given MAG ID,
                    or if all contigs are below the length cutoff.
    """

    logging.info(f"Reading contig lengths from: {contigs_file}")
    length_df = pd.read_csv(contigs_file, sep="\t")

    if "contig" not in length_df.columns or "length" not in length_df.columns:
        logging.error("contig_lengths file must have columns: contig, length")
        sys.exit(1)

    # Filter to the selected MAG
    mag_subset_df = length_df[length_df["MAG"] == mag_id]

    if mag_subset_df.empty:
        logging.error(f"No contigs found for MAG {mag_id}. Exiting.")
        sys.exit(42)

    logging.info(f"Found {len(mag_subset_df):,} contigs belonging to MAG {mag_id}.")

    # Filter by length_cutoff
    mag_subset_df = mag_subset_df[mag_subset_df["length"] >= length_cutoff]
    logging.info(
        f"After filtering contigs < {length_cutoff:,} kbp, we have {len(mag_subset_df)} contigs left for MAG {mag_id}."
    )

    if mag_subset_df.empty:
        logging.warning("All contigs for this MAG are below length cutoff. Exiting.")
        sys.exit(42)

    valid_contigs = set(mag_subset_df["contig"])

    return valid_contigs


def filter_significance(significance_file, valid_contigs):
    """
    Filters a significance table to include only rows with contigs present in the valid_contigs list.

    Parameters:
        significance_file (str): Path to the significance table file in tab-separated values format.
        valid_contigs (list): List of contigs that are considered valid.

    Returns:
        pandas.DataFrame: Filtered significance table containing only rows with valid contigs.

    Logs:
        Logs the number of rows before and after filtering.
        Logs a warning and exits with status code 42 if no rows remain after filtering.
    """

    logging.info(f"Reading significance table: {significance_file}")
    sig_df = pd.read_csv(significance_file, sep="\t")

    original_count = len(sig_df)
    # Keep only the contigs that passed the length filter
    sig_df = sig_df[sig_df["contig"].isin(valid_contigs)]
    filtered_count = len(sig_df)

    logging.info(
        f"Significance table: filtered from {original_count:,} -> {filtered_count:,} rows "
        f"based on valid contigs."
    )

    if sig_df.empty:
        logging.warning(
            "No rows remain in significance table after filtering. Exiting."
        )
        sys.exit(42)

    return sig_df


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Perform a sliding window significance analysis for a specific MAG, filtering out short contigs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--mag_id",
        required=True,
        type=str,
        help="MAG ID of interest (e.g., 'SLG1007_DASTool_bins_35').",
    )

    parser.add_argument(
        "--contig_lengths",
        required=True,
        help="Path to a table of contigs across MAGs (col: contig, length).",
    )

    parser.add_argument(
        "--significance_table",
        required=True,
        type=str,
        help="Path to the significance data (TSV/CSV) that has contig, position, p_value_* columns.",
    )

    parser.add_argument(
        "--length_cutoff",
        type=int,
        default=10,
        help="Minimum contig length (in kbp) for analysis.",
    )

    parser.add_argument(
        "--window_size", type=int, default=100, help="Window size in bp."
    )

    parser.add_argument(
        "--p_value_threshold",
        type=float,
        default=0.05,
        help="p-value threshold for significance.",
    )

    parser.add_argument(
        "--output_table",
        required=True,
        type=str,
        help="Path to output the final sliding-window table (TSV).",
    )

    args = parser.parse_args()

    start_time = time.time()
    # Create a set of contigs that pass the length filter
    valid_contigs = filter_contigs(
        args.contig_lengths, args.mag_id, args.length_cutoff * 1000
    )

    sig_df = filter_significance(args.significance_table, valid_contigs)

    # Extract test columns
    logging.info("Extracting test columns...")
    test_columns_dict = extract_test_columns(sig_df)

    logging.info(
        f"Computing sliding window significance with window_size={args.window_size}, "
        f"p_value_threshold={args.p_value_threshold}"
    )
    result_df = sliding_window_significance(
        sig_df,
        test_columns_dict=test_columns_dict,
        window_size=args.window_size,
        p_value_threshold=args.p_value_threshold,
    )

    logging.info(f"Sliding window result has {len(result_df):,} rows.")

    # Write final table
    logging.info(f"Writing final table to {args.output_table}")
    result_df.to_csv(args.output_table, sep="\t", index=False)

    end_time = time.time()
    logging.info(f"Done! Time taken {end_time - start_time:.2f} seconds.")


if __name__ == "__main__":
    main()
