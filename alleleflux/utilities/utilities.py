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


def calculate_mag_sizes(fasta_file, mag_mapping_file):
    """
    Calculate the sizes of Metagenome-Assembled Genomes (MAGs) from a FASTA file.

    This function parses a given FASTA file and calculates the total size of each MAG
    by summing the lengths of its contigs. The MAG ID is determined either from a mapping file
    or by extracting from the contig ID (default: splits by '.fa').

    Parameters:
        fasta_file (str): Path to the input FASTA file containing contig sequences.
        mag_mapping_file (str): Path to a file mapping contigs to MAG IDs.

    Returns:
        dict: A dictionary where keys are MAG IDs and values are the total sizes of the MAGs.
    """
    logging.info("Parsing FASTA file to calculate MAG size.")

    # Load MAG mapping if provided
    mag_mapping = load_mag_mapping(mag_mapping_file)

    mag_sizes = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_id = record.id
        # Extract MAG ID from contig ID using mapping or default method
        mag_id = extract_mag_id(contig_id, mag_mapping)
        # Accumulate the length of the contig to the MAG's total size
        mag_sizes[mag_id] += len(record.seq)

    logging.info(f"Calculated sizes for {len(mag_sizes)} MAGs")
    return mag_sizes


def load_mag_mapping(mapping_file):
    """
    Load a MAG-to-contig mapping file.

    Parameters:
        mapping_file (str): Path to the mapping file (tab-separated values with columns: mag_id, contig_id).

    Returns:
        dict: A dictionary where keys are contig names and values are MAG IDs.

    Raises:
        ValueError: If the mapping file is not provided, missing required columns, or cannot be read.
    """
    if not mapping_file:
        raise ValueError(
            "MAG mapping file is required. Please provide a valid path to a MAG-to-contig mapping file."
        )

    logging.info(f"Loading MAG-to-contig mapping from: {mapping_file}")
    df = pd.read_csv(mapping_file, sep="\t")

    required_columns = {"mag_id", "contig_id"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(
            f"Missing columns in MAG mapping file: {missing_columns}. "
            f"Required columns are: {required_columns}"
        )

    # Create a dictionary mapping contig names to MAG IDs
    mapping = dict(zip(df["contig_id"], df["mag_id"]))
    logging.info(f"Loaded mapping for {len(mapping):,} contigs")
    return mapping


def extract_mag_id(contig_id, mag_mapping):
    """
    Extract MAG ID from a contig ID using the provided mapping.

    Parameters:
        contig_id (str): The contig ID from which to extract the MAG ID.
        mag_mapping (dict): Dictionary mapping contig IDs to MAG IDs.

    Returns:
        str: The mapped MAG ID.

    Raises:
        KeyError: If the contig ID is not found in the mapping.
        TypeError: If mag_mapping is None.
    """
    if mag_mapping is None:
        raise TypeError(
            "mag_mapping cannot be None. Check if the mapping was properly initialized."
        )

    if contig_id not in mag_mapping:
        raise KeyError(
            f"Contig ID '{contig_id}' not found in the MAG mapping. Please ensure your mapping file contains all contigs."
        )

    return mag_mapping[contig_id]


def load_mag_metadata_file(
    mag_metadata_file, mag_id, breath_threshold, data_type="longitudinal"
):
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

    if data_type == "longitudinal":
        required_columns = {
            "sample_id",
            "file_path",
            "subjectID",
            "group",
            "time",
            "replicate",
        }
    elif data_type == "single":  # single data
        required_columns = {"sample_id", "file_path", "subjectID", "group", "replicate"}

    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing columns in mag metadata file: {missing_columns}")

    # Convert columns to string (only for the ones we require)
    common_cast = {
        "sample_id": str,
        "file_path": str,
        "subjectID": str,
        "group": str,
        "replicate": str,
    }
    if data_type == "longitudinal":
        cast_dict = {**common_cast, "time": str}
        df = df.astype(cast_dict)
        metadata_dict = df.set_index("sample_id")[
            ["group", "time", "subjectID", "replicate"]
        ].to_dict(orient="index")
    else:
        df = df.astype(common_cast)
        metadata_dict = df.set_index("sample_id")[
            ["group", "subjectID", "replicate"]
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


def extract_relevant_columns(df, capture_str, lmm_format=False):
    """
    Extract columns from DataFrame based on a capture string pattern.

    Parameters:
        df (pandas.DataFrame): DataFrame containing columns to extract.
        capture_str (str): String pattern to look for in column names.
        lmm_format (bool, optional): If True, handle LMM-style column names (e.g. 'A_p_value')
                                     instead of the default format ('A_frequency_p_value_Wilcoxon').

    Returns:
        dict: Dictionary mapping test names to lists of relevant column names.
    """
    if lmm_format:
        # Handle LMM.py output format (e.g., 'A_p_value', 'T_p_value', etc.)
        nucleotides = ["A", "T", "G", "C"]
        str_columns = [
            col
            for col in df.columns
            if col.endswith("_p_value") and col.split("_")[0] in nucleotides
        ]
        if not str_columns:
            raise ValueError("No LMM p-value columns found.")

        test_columns_dict = {"lmm": str_columns}
        logging.info(f"Detected LMM p-value columns: {str_columns}")
    else:
        test_columns_dict = {}
        str_columns = [col for col in df.columns if capture_str in col]
        if not str_columns:
            raise ValueError(f"No columns found containing '{capture_str}'")

        for col in str_columns:
            # Everything after 'capture_str' is part of the test name
            test_name = col.split(capture_str)[-1]

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

        # Check for NaNs in p-value columns and log them instead of raising an error
        nan_counts = subdf[test_cols].isnull().sum().sum()
        if nan_counts > 0:
            logging.warning(
                f"Found {nan_counts} NaN values in {test_name} p-value columns. NAs will be skipped in calculations."
            )

        # The `dropna=False` ensures that groups with NaN values are included as a separate group
        grouped = subdf.groupby(group_by_column, dropna=False)

        # Total sites = number of present sites per group
        total_sites_per_group = grouped.size()

        # Identify significant sites
        significant_col = f"is_significant_{test_name}"
        # This returns False for cells with NAs
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


def read_gtdb(gtdb_fpath):
    gtdb_df = pd.read_csv(
        gtdb_fpath,
        sep="\t",
        usecols=["user_genome", "classification"],
    )
    gtdb_df = gtdb_df.rename(columns={"user_genome": "MAG_ID"})
    taxon_data = gtdb_df["classification"].apply(parse_classification)
    taxon_df = pd.DataFrame(taxon_data.tolist())
    gtdb_df = pd.concat([gtdb_df, taxon_df], axis=1)
    gtdb_df = gtdb_df.drop(columns=["classification"])
    return gtdb_df


def parse_classification(classification_str):
    """
    Parses a taxonomic classification string and returns a dictionary with taxonomic ranks.

    Parameters:
        classification_str (str): A semicolon-separated string containing taxonomic classifications
                                  with prefixes (e.g., "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria").

    Returns:
        dict: A dictionary where keys are taxonomic ranks (e.g., "domain", "phylum", "class", etc.)
              and values are the corresponding names from the classification string. If a rank is
              not specified in the input string, its value will be "unclassified".
    """
    # Define taxonomic ranks and prefixes
    taxonomic_ranks = [
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    rank_prefixes = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
    # Set all taxonomic ranks to unclassified
    taxon_dict = {rank: "unclassified" for rank in taxonomic_ranks}
    taxa = classification_str.split(";")
    for taxon in taxa:
        for prefix, rank in zip(rank_prefixes, taxonomic_ranks):
            if taxon.startswith(prefix):
                # remove the prefix and get the name
                name = taxon.replace(prefix, "").strip()
                if name == "":
                    name = "unclassified"
                taxon_dict[rank] = name
                # Exit the loop to avoid further unnecessary iterations
                break
    return taxon_dict
