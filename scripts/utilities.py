#!/usr/bin/env python
import logging
from collections import defaultdict
from functools import reduce

import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)


def calculate_mag_sizes(fasta_file):
    """
    Calculate the sizes of Metagenome-Assembled Genomes (MAGs) from a FASTA file.

    This function parses a given FASTA file and calculates the total size of each MAG
    by summing the lengths of its contigs. The MAG ID is extracted from the contig ID,
    which is assumed to be the part of the contig ID before '.fa'.

    Parameters:
        fasta_file (str): Path to the input FASTA file containing contig sequences.

    Returns:
        dict: A dictionary where keys are MAG IDs and values are the total sizes of the MAGs.
    """
    logging.info("Parsing FASTA file to calculate MAG size.")
    mag_sizes = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_id = record.id
        # Extract MAG ID from contig ID
        # The MAG ID is everything before '.fa' in the contig ID
        mag_id = contig_id.split(".fa")[0]
        # Accumulate the length of the contig to the MAG's total size
        mag_sizes[mag_id] += len(record.seq)
    return mag_sizes


def load_mag_metadata_file(mag_metadata_file, mag_id, breath_threshold):
    """
    Load MAG metadata from a file and process it.

    Parameters:
        mag_metadata_file (str): Path to the MAG metadata file (tab-separated values).
        mag_id (str): The MAG identifier to be associated with each sample.
        breath_threshold (float): The breath threshold value to be associated with each sample.

    Returns:
        tuple: A tuple containing:
            - metadata_dict (dict): A dictionary where keys are sample IDs and values are dictionaries with keys 'group' and 'subjectID'.
            - sample_files_with_mag_id (list): A list of tuples, each containing (sample_id, file_path, mag_id, breath_threshold).

    Raises:
        ValueError: If the MAG metadata file is missing required columns.
    """

    logging.info(f"Loading MAG metadata from file: {mag_metadata_file}")
    df = pd.read_csv(mag_metadata_file, sep="\t")

    # Ensure required columns are present
    required_columns = {
        "sample_id",
        "file_path",
        "subjectID",
        "group",
        "time",
        "replicate",
    }
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing columns in mag metadata file: {missing_columns}")

    # Convert columns to string
    df = df.astype(
        {
            "sample_id": str,
            "file_path": str,
            "subjectID": str,
            "group": str,
            "time": str,
            "replicate": str,
        }
    )

    # Build metadata_dict
    metadata_dict = df.set_index("sample_id")[
        ["group", "time", "subjectID", "replicate"]
    ].to_dict(orient="index")

    # Build sample_files_with_mag_id: list of (sample_id, file_path, mag_id, breath_threshold)
    sample_files_with_mag_id = (
        df[["sample_id", "file_path"]]
        .apply(
            lambda row: (row["sample_id"], row["file_path"], mag_id, breath_threshold),
            axis=1,
        )
        .tolist()
    )

    return metadata_dict, sample_files_with_mag_id


def extract_relevant_columns(df, capture_str):
    str_columns = [col for col in df.columns if capture_str in col]

    test_columns_dict = {}

    for col in str_columns:
        if capture_str in col:
            # Everything after 'capture_str' is part of the test name
            test_name = col.split(capture_str)[-1]
        else:
            raise ValueError(
                f"Column {col} does not contain {capture_str} in expected format."
            )

        test_columns_dict.setdefault(test_name, []).append(col)
    logging.info(f"Detected tests: {list(test_columns_dict.keys())}")
    return test_columns_dict


def calculate_score(df, test_columns_dict, group_by_column, p_value_threshold=0.05):

    results = []

    for test_name, test_cols in test_columns_dict.items():
        # Create a subset DataFrame with only the grouping column and test columns
        subset_cols = [group_by_column] + test_cols
        subdf = df[subset_cols].copy()
        subdf.dropna(subset=test_cols, how="all", inplace=True)

        # Check for NaNs in all p-value columns and drop
        if subdf[test_cols].isnull().values.any():
            raise ValueError("NaNs found in p-value column.")

        # The `dropna=False` ensures that groups with NaN values are included as a separate group
        grouped = subdf.groupby(group_by_column, dropna=False)

        # Total sites = number of present sites per group
        total_sites_per_group = grouped.size()

        # Identify significant sites
        significant_col = f"is_significant_{test_name}"
        subdf[significant_col] = subdf[test_cols].lt(p_value_threshold).any(axis=1)

        # Number of significant sites per group
        significant_sites_per_group = grouped[significant_col].sum()

        # Percentage of significant sites per group
        percentage_significant = (
            significant_sites_per_group / total_sites_per_group
        ) * 100

        test_result = pd.DataFrame(
            {
                group_by_column: total_sites_per_group.index,
                f"total_sites_per_group_{test_name}": total_sites_per_group.values,
                f"significant_sites_per_group_{test_name}": significant_sites_per_group.values,
                f"score_{test_name} (%)": percentage_significant.values,
            }
        )

        results.append(test_result)

    # Merge all test results into a single DataFrame
    final_results = reduce(
        lambda left, right: pd.merge(left, right, on=group_by_column, how="outer"),
        results,
    )
    return final_results
