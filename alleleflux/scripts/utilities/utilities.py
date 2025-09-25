#!/usr/bin/env python
import logging
from collections import defaultdict
from functools import reduce

import pandas as pd
from Bio import SeqIO

# Set up logger for this module
logger = logging.getLogger(__name__)


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
    logger.info("Parsing FASTA file to calculate MAG size.")

    # Load MAG mapping if provided
    mag_mapping = load_mag_mapping(mag_mapping_file)

    mag_sizes = defaultdict(int)
    for record in SeqIO.parse(fasta_file, "fasta"):
        contig_id = record.id
        # Extract MAG ID from contig ID using mapping or default method
        mag_id = extract_mag_id(contig_id, mag_mapping)
        # Accumulate the length of the contig to the MAG's total size
        mag_sizes[mag_id] += len(record.seq)

    logger.info(f"Calculated sizes for {len(mag_sizes)} MAGs")
    return mag_sizes


def build_contig_length_index(fasta_file, mag_mapping_file):
    """
    Build an index of contig lengths limited to contigs present in the mapping file.

    Parameters
    ----------
    fasta_file : str
        Path to FASTA with contig sequences.
    mag_mapping_file : str
        Tab-separated mapping file with columns: mag_id, contig_id.

    Returns
    -------
    dict
        contig_id -> length (only for contigs referenced in mapping file)
    """
    logger.info("Building contig length index...")
    mag_mapping_df = pd.read_csv(mag_mapping_file, sep="\t")
    required_columns = {"mag_id", "contig_id"}
    missing = required_columns - set(mag_mapping_df.columns)
    if missing:
        raise ValueError(
            f"Missing columns in MAG mapping file: {missing}. Required: {required_columns}"
        )

    contig_whitelist = set(mag_mapping_df["contig_id"].tolist())
    contig_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in contig_whitelist:
            contig_lengths[record.id] = len(record.seq)

    logger.info(f"Indexed {len(contig_lengths):,} contigs with lengths for coverage weighting.")
    return contig_lengths


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

    logger.info(f"Loading MAG-to-contig mapping from: {mapping_file}")
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
    logger.info(f"Loaded mapping for {len(mapping):,} contigs")
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

    logger.info(f"Loading MAG metadata from file: {mag_metadata_file}")
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


def extract_relevant_columns(df, capture_str):
    """
    Extract columns from DataFrame based on a capture string pattern.

    Parameters:
        df (pandas.DataFrame): DataFrame containing columns to extract.
        capture_str (str): String pattern to look for in column names.

    Returns:
        dict: Dictionary mapping test names to lists of relevant column names.
    """

    test_columns_dict = {}
    str_columns = [col for col in df.columns if capture_str in col]
    if not str_columns:
        raise ValueError(f"No columns found containing '{capture_str}'")

    for col in str_columns:
        # Everything after 'capture_str' is part of the test name
        test_name = col.split(capture_str)[-1]

        test_columns_dict.setdefault(test_name, []).append(col)
    logger.info(f"Detected tests: {list(test_columns_dict.keys())}")
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
            logger.warning(
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


def load_and_filter_data(
    input_df_path: str,
    preprocessed_df_path: str,
    mag_id: str,
    dtype_map: dict,
    group_to_analyze: str = None,
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
    # Read the preprocessed positions for filtering
    logger.info(f"Loading preprocessed positions from {preprocessed_df_path}")
    preprocessed_df = pd.read_csv(
        preprocessed_df_path,
        sep="\t",
        usecols=["contig", "position", "group"],
    )
    logger.info(f"Loaded {preprocessed_df.shape[0]:,} preprocessed positions.")

    # Verify the requested group exists in the data. Required for LMM parallelism score.
    if group_to_analyze and group_to_analyze not in preprocessed_df["group"].unique():
        raise ValueError(f"Group '{group_to_analyze}' not found in the data")

    preprocessed_df.drop(columns=["group"], inplace=True)

    # Read the raw allele count data
    logger.info(f"Loading raw count data from {input_df_path}")
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
    logger.info(f"Loaded {raw_counts_df.shape[0]:,} rows of raw count data.")
    # Detailed stats are debug level for performance
    logger.debug(
        f"Distinct contigs: {raw_counts_df.index.get_level_values('contig').nunique():,}",
    )
    logger.debug(
        f"Distinct positions: {raw_counts_df.index.nunique():,}",
    )
    # Get unique (contig, position) pairs to keep
    valid_positions = preprocessed_df.drop_duplicates(subset=["contig", "position"])
    valid_positions_index = pd.MultiIndex.from_frame(valid_positions)

    # Filter Raw Counts by Preprocessed Positions
    logger.info("Filtering raw counts to include only preprocessed positions.")

    # Keep only rows whose (contig, position) index is in the valid_positions_index
    filtered_counts_df = raw_counts_df[
        raw_counts_df.index.isin(valid_positions_index)
    ].reset_index()

    logger.info(
        f"Filtered to {filtered_counts_df.shape[0]:,} rows after applying preprocessed positions."
    )

    if filtered_counts_df.empty:
        raise ValueError(
            f"No positions remaining after filtering by preprocessed data for MAG {mag_id}. Exiting."
        )

    logger.debug(
        f"{len(filtered_counts_df['contig'].unique()):,} contigs and "
        f"{len(filtered_counts_df.groupby(['contig', 'position'], dropna=False)):,} unique positions remaining after filtering."
    )
    return filtered_counts_df
