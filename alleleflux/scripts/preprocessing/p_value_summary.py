#!/usr/bin/env python3
"""
Compile and Filter Significant Sites from AlleleFlux Test Results
=================================================================

This script processes and consolidates the output from various AlleleFlux
significance tests. It automatically discovers result files within a given
directory, includes all sites in the analysis, and applies FDR-BH correction
to control for multiple testing.

The script intelligently handles different p-value column naming conventions
from various tests (e.g., LMM, CMH, t-tests) and applies FDR-BH correction
to control for multiple testing.

Key Features:
- Automatic discovery of test result files from a single input directory.
- Intelligent identification of p-value columns for different test types.
- Inclusion of all sites for comprehensive FDR correction analysis.
- Consolidation of results from multiple tests and MAGs into one table.
- FDR-BH correction applied to each test type's concatenated results.
- Optional grouping by MAG ID for more stringent FDR correction.
- Optional generation of a human-readable summary report.

Usage:
    python compile_test_results.py --input-dir /path/to/alleleflux_results
                                 --outdir /path/to/output
                                 --timepoints pre_post
                                 --test-types lmm cmh paired_tTest
                                 --group-by-mag-id
"""

import argparse
import logging
import sys
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

from alleleflux.scripts.utilities.utilities import extract_relevant_columns

# Set up logger for this script
logger = logging.getLogger(__name__)
# LOGGER: logging.Logger = logging.getLogger(__name__)


@contextmanager
def suppress_logger(logger_name, level=logging.WARNING):
    """Context manager to temporarily suppress a logger's output."""
    target_logger = logging.getLogger(logger_name)
    original_level = target_logger.level
    target_logger.setLevel(level)
    try:
        yield
    finally:
        target_logger.setLevel(original_level)


def find_test_files(
    input_dir: Path, test_types: list[str], timepoint_label: str
) -> dict[str, list[Path]]:
    """
    Finds all AlleleFlux test result files in the input directory for the specified test types and timepoint label.

    Args:
        input_dir: The directory containing AlleleFlux result folders.
        test_types: A list of test types to search for (e.g., ['lmm', 'cmh']).
        timepoint_label: The timepoint label to include (e.g., 'pre_post' or 'pre').

    Returns:
        A dictionary mapping each test type to a list of its result file paths.
    """
    test_files = {test_type: [] for test_type in test_types}
    logger.info(f"Searching for test results in: {input_dir}")

    for test_type in test_types:
        # Search for directories that start with the test type name
        # e.g., 'lmm', 'two_sample_paired'
        for test_dir in input_dir.glob(f"{test_type}_{timepoint_label}-*"):
            if test_dir.is_dir():
                found_files = list(test_dir.glob("*.tsv.gz"))

                if found_files:
                    test_files[test_type].extend(found_files)
                    logger.info(
                        f"  Found {len(found_files)} file(s) for test '{test_type}' in {test_dir.name}"
                    )

    return test_files


def process_results_file(
    filepath: Path, test_type: str, timepoint_label: str
) -> dict[str, pd.DataFrame]:
    """
    Processes a single result file, consolidating all sub-tests into a single table.
    The output filename is determined by the test type and timepoint label, not the
    sub-test name. A new 'test_type' column is added to the DataFrame to specify
    the sub-test, and a 'group_analyzed' column is added for relevant tests.
    Args:
        filepath: Path to the AlleleFlux result file.
        test_type: The name of the test being processed.
        timepoint_label: The timepoint label to include.
    Returns:
        A dictionary mapping a generated file group key to a consolidated DataFrame
        for the given file, or an empty dictionary if no sites are found.
    """
    # Load the data
    df = pd.read_csv(filepath, sep="\t", compression="gzip", low_memory=False)
    if df.empty:
        logger.warning(f"File {filepath} for test '{test_type}' is empty. Skipping.")
        return {}
    # Identify all possible p-value columns, grouped by the sub-test
    p_value_cols_by_test = extract_relevant_columns(df, capture_str="p_value_")
    if not p_value_cols_by_test:
        logger.warning(f"No p-value columns found in {filepath}. Skipping.")
        return {}
    # Extract MAG ID once for the entire file
    if "mag_id" not in df.columns:
        mag_id = extract_mag_id_from_filepath(filepath, test_type)
    else:
        # Check if all rows have the same MAG ID
        unique_mag_ids = df["mag_id"].unique()
        if len(unique_mag_ids) != 1:
            raise ValueError(
                f"Multiple MAG IDs found in file {filepath}: {unique_mag_ids}. "
                f"Expected exactly one MAG ID per file."
            )
        mag_id = unique_mag_ids[0]
    # Add common metadata columns once
    df["mag_id"] = mag_id
    df["source_file"] = filepath.name
    df["period"] = timepoint_label
    all_subtest_dfs = []
    # Process each sub-test
    for sub_test_name, p_value_cols in p_value_cols_by_test.items():
        # Work on a copy to avoid modifying the original DataFrame in the loop
        sub_df = df.copy()
        # Calculate minimum p-value for this sub-test
        sub_df["min_p_value"] = sub_df[p_value_cols].min(axis=1)
        # Set test_type and group_analyzed based on test_type and sub_test_name
        new_test_type = f"{test_type}_{sub_test_name}"  # Default
        if test_type in ["lmm", "lmm_across_time"]:
            new_test_type = "LMM"
            if test_type == "lmm_across_time" and "group_analyzed" in df.columns:
                sub_df["group_analyzed"] = df["group_analyzed"]
        elif test_type in ["cmh", "cmh_across_time"]:
            new_test_type = "CMH"
            if test_type == "cmh_across_time" and "analyzed_group" in df.columns:
                sub_df["group_analyzed"] = df["analyzed_group"]
        elif test_type == "single_sample":
            # Handles "tTest_{group}" and "Wilcoxon_{group}"
            if sub_test_name.startswith("tTest_"):
                new_test_type = "single_sample_tTest"
                sub_df["group_analyzed"] = sub_test_name[len("tTest_") :]
            elif sub_test_name.startswith("Wilcoxon_"):
                new_test_type = "single_sample_Wilcoxon"
                sub_df["group_analyzed"] = sub_test_name[len("Wilcoxon_") :]
            else:
                raise ValueError(
                    f"Unexpected sub-test name '{sub_test_name}' for single_sample test type in {filepath}. "
                )

        sub_df["test_type"] = new_test_type
        # if "group_analyzed" not in sub_df:
        #     sub_df["group_analyzed"] = None

        # Select final columns
        output_cols = [
            "period",
            "mag_id",
            "contig",
            "position",
            "gene_id",
            "test_type",
            "group_analyzed",
            "min_p_value",
            "source_file",
        ]
        # Store results (only columns that exist) for this sub-test
        available_cols = [col for col in output_cols if col in sub_df.columns]
        all_subtest_dfs.append(sub_df[available_cols])

    if not all_subtest_dfs:
        logger.warning(f"No sites found in {filepath} for test '{test_type}'.")
        return {}

    # Create a filename-friendly key that groups results by test type and timepoint
    combined_key = f"{test_type}_{timepoint_label}"
    # Consolidate all sub-test data from this file into one DataFrame
    consolidated_df = pd.concat(all_subtest_dfs, ignore_index=True)

    return {combined_key: consolidated_df}


def extract_mag_id_from_filepath(filepath: Path, test_type: str) -> str:
    """
    Extract MAG ID from the filename using AlleleFlux naming conventions.

    Args:
        filepath: Path to the result file
        test_type: The test type being processed (e.g., 'lmm_across_time', 'two_sample_paired', 'cmh_across_time')

    Returns:
        MAG ID extracted from filename
    """
    filename = filepath.name

    # Look for the test_type pattern in the filename
    if f"_{test_type}" in filename:
        # Split on the test_type pattern and take the first part as MAG ID
        mag_id = filename.split(f"_{test_type}")[0]
        return mag_id


def main():
    """Main function to drive the script."""
    parser = argparse.ArgumentParser(
        description="Compile and filter significant sites from AlleleFlux test results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing subdirectories of AlleleFlux test results.",
    )
    parser.add_argument(
        "--p-value-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for reference (currently unused as all sites are included for FDR correction).",
    )
    parser.add_argument(
        "--timepoints",
        type=str,
        required=True,
        help="A single timepoint combination label to include (e.g., 'pre_post' or 'pre').",
    )
    parser.add_argument(
        "--test-types",
        nargs="+",
        default=[
            "two_sample_unpaired",
            "two_sample_paired",
            "lmm",
            "lmm_across_time",
            "cmh",
            "cmh_across_time",
            "single_sample",
        ],
        help="Space-separated list of test types to include (e.g., lmm cmh paired_tTest) (default: all).",
    )
    parser.add_argument(
        "--outdir", required=True, type=Path, help="Directory to save the output files."
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="p_value_summary",
        help="Prefix for the output file names.",
    )
    parser.add_argument(
        "--fdr-group-by-mag-id",
        action="store_true",
        help="Apply FDR correction separately for each MAG ID in addition to test_type and group_analyzed.",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # --- Main Workflow ---
    if not args.input_dir.exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    # 1. Discover all relevant result files
    test_files = find_test_files(args.input_dir, args.test_types, args.timepoints)
    total_files = sum(len(files) for files in test_files.values())
    if total_files == 0:
        logger.error(
            "No result files found for the specified test types and timepoints. Exiting."
        )
        sys.exit(1)
    logger.info(f"Found a total of {total_files} files to process.")

    # 2. Process each file and collect significant sites by sub-test
    all_significant_sites_by_subtest = {}
    for test_type, file_paths in test_files.items():
        if not file_paths:
            logger.warning(f"No files found for test type '{test_type}'. Skipping.")
            continue

        logger.info(f"Processing {len(file_paths)} files for test type: '{test_type}'")
        for filepath in tqdm(
            file_paths,
            desc=f"Processing {test_type}",
            unit="file",
            total=len(file_paths),
        ):
            # Suppress utilities logger to avoid interrupting tqdm
            with suppress_logger("alleleflux.scripts.utilities.utilities"):
                significant_sites_by_subtest = process_results_file(
                    filepath, test_type, args.timepoints
                )

            # Organize results by sub-test with combined key
            for combined_key, significant_df in significant_sites_by_subtest.items():
                # Initialize the list for this combined key if it doesn't exist,
                # this prevents KeyError, when initializing the list for the first time for a combined key
                if combined_key not in all_significant_sites_by_subtest:
                    all_significant_sites_by_subtest[combined_key] = []
                all_significant_sites_by_subtest[combined_key].append(significant_df)

    # 3. Consolidate, sort, and save separate tables for each sub-test
    # Process each sub-test separately
    for combined_key, subtest_dataframes in all_significant_sites_by_subtest.items():
        logger.info(f"Consolidating significant sites for: {combined_key}")

        # Concatenate all dataframes for this sub-test
        if not subtest_dataframes:
            logger.warning(f"No dataframes for key '{combined_key}', skipping.")
            continue
        final_df = pd.concat(subtest_dataframes, ignore_index=True)

        # Define grouping columns for FDR correction. Correction is applied per group.
        grouping_cols = ["test_type"]
        if "group_analyzed" in final_df.columns:
            grouping_cols.append("group_analyzed")

        # Group by mag_id is needed.
        if args.fdr_group_by_mag_id:
            grouping_cols.append("mag_id")

        logger.info(
            f"Applying FDR-BH correction grouping cols: {grouping_cols} within {combined_key}..."
        )

        # Apply FDR-BH correction within each group using transform.
        final_df["q_value"] = final_df.groupby(grouping_cols)["min_p_value"].transform(
            lambda x: multipletests(x, method="fdr_bh", alpha=args.p_value_threshold)[1]
        )

        # Create output filename using the combined key
        subtest_output_file = args.outdir / f"{args.prefix}_{combined_key}.tsv"

        args.outdir.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(subtest_output_file, sep="\t", index=False)
        logger.info(
            f"Successfully wrote {len(final_df):,} sites for {combined_key} to: {subtest_output_file}"
        )

    logger.info("Analysis complete.")


if __name__ == "__main__":
    main()
