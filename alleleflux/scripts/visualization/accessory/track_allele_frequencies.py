#!/usr/bin/env python3
"""
Track allele frequencies for a specified MAG based on "anchor" alleles.

This script takes an "anchor file" (e.g., the output from the
terminal_nucleotide_analysis.py script) which defines the allele of interest
for each genomic position.

It then scans the profile files for *all* samples (as defined in the
metadata) and calculates the frequency of that specific anchor allele
in every sample.

This produces a "wide-form" table (sites x samples) and a long-form table that is "tidy"
and ready for analysis in R, with one row per site-allele combination.
"""

import argparse
import logging
import multiprocessing
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# Global constant for nucleotide symbols
NUCLEOTIDES = ["A", "C", "G", "T"]


def validate_inputs(args):
    """Validate input arguments and create output directory."""
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.anchor_file.exists():
        logger.error(f"Anchor file not found: {args.anchor_file}")
        raise FileNotFoundError(f"Anchor file not found: {args.anchor_file}")
    if not args.metadata.exists():
        logger.error(f"Metadata file not found: {args.metadata}")
        raise FileNotFoundError(f"Metadata file not found: {args.metadata}")


def load_and_validate_metadata(metadata_path):
    """Load and validate metadata file with required columns."""
    logger.info(f"Loading metadata from: {metadata_path}")
    metadata_df = pd.read_csv(metadata_path, sep="\t")

    # Validate required metadata columns
    required_meta_cols = [
        "sample_id",
        "sample_profile_dir",
        "group",
        "time",
        "subjectID",
    ]
    missing_cols = set(required_meta_cols) - set(metadata_df.columns)
    if missing_cols:
        logger.error(f"Metadata is missing required columns: {missing_cols}")
        raise ValueError(f"Metadata is missing required columns: {missing_cols}")

    # Log optional columns
    optional_cols = ["day", "replicate"]
    present_optional = [col for col in optional_cols if col in metadata_df.columns]
    if present_optional:
        logger.info(f"Optional columns present: {present_optional}")
    missing_optional = [col for col in optional_cols if col not in metadata_df.columns]
    if missing_optional:
        logger.info(f"Optional columns not present: {missing_optional}")

    return metadata_df


def load_and_prepare_anchor_sites(anchor_file_path, mag_id, anchor_column):
    """Load anchor file and prepare anchor alleles with explode logic."""
    logger.info(f"Loading anchor alleles for {mag_id} from: {anchor_file_path}")
    anchor_sites_df = pd.read_csv(anchor_file_path, sep="\t")

    # Filter for the specific MAG we are processing
    anchor_mag_df = anchor_sites_df[anchor_sites_df["mag_id"] == mag_id].copy()

    if anchor_mag_df.empty:
        logger.error(f"No sites found for MAG '{mag_id}' in anchor file.")
        raise ValueError(f"No sites found for MAG '{mag_id}' in anchor file.")

    if anchor_column not in anchor_mag_df.columns:
        logger.error(f"Anchor column '{anchor_column}' not found in anchor file.")
        raise ValueError(f"Anchor column '{anchor_column}' not found.")

    # Count NaN values in anchor_column before processing
    nan_count = anchor_mag_df[anchor_column].isna().sum()
    if nan_count > 0:
        logger.info(
            f"Found {nan_count} sites with NaN values in '{anchor_column}' for MAG '{mag_id}'"
        )

    # Handle NaNs before splitting, cast to string
    anchor_mag_df[anchor_column] = anchor_mag_df[anchor_column].fillna("").astype(str)

    # Create the new 'anchor_allele' column by splitting the anchor_column
    anchor_mag_df["anchor_allele"] = anchor_mag_df[anchor_column].str.split(",")

    # Explode the DataFrame, one row for each allele in the list
    exploded_df = anchor_mag_df.explode("anchor_allele")

    # Filter out empty strings that came from NaNs
    exploded_df = exploded_df[exploded_df["anchor_allele"] != ""].copy()
    # Filter anchors to A/C/G/T and warn, so you don’t carry junk values into workers
    exploded_df["anchor_allele"] = exploded_df["anchor_allele"].str.strip().str.upper()
    bad = ~exploded_df["anchor_allele"].isin(NUCLEOTIDES)
    if bad.any():
        logger.warning(f"Dropping {int(bad.sum())} anchors not in {NUCLEOTIDES}")
        exploded_df = exploded_df[~bad].copy()

    if exploded_df.empty:
        logger.error(f"No valid anchor alleles found for MAG '{mag_id}'")
        raise ValueError(f"No valid anchor alleles found for MAG '{mag_id}'")

    # Set index for easy merging later
    key_cols = ["contig", "position", "anchor_allele"]

    # Check for duplicate keys
    if exploded_df.duplicated(subset=key_cols).any():
        raise ValueError(
            f"Found duplicate (contig, position, allele) entries for {mag_id}. "
            "This can happen if the anchor file has errors. Using first entry."
        )

    final_anchor_info_df = exploded_df.set_index(key_cols)

    # Create the anchor_dict for the worker
    anchor_dict = (
        final_anchor_info_df.reset_index()
        .groupby(["contig", "position"])["anchor_allele"]
        .apply(list)
        .to_dict()
    )

    if not anchor_dict:
        logger.error("Anchor dictionary is empty. No sites to process.")
        raise ValueError("No valid anchor sites to process.")

    logger.info(f"Loaded {len(anchor_dict):,} sites to check for {mag_id}.")
    return final_anchor_info_df, anchor_dict


def process_single_sample(sample_data, anchor_dict, mag_id, min_cov_per_site):
    """Worker function to process a single sample's profile file.

    For a given sample, this function reads its profile file and calculates
    the frequency of all specified anchor alleles.

    Parameters:
    -----------
    sample_data : tuple
        A tuple containing (sample_id, sample_profile_dir) from the metadata.
    anchor_dict : dict
        A dictionary mapping (contig, position) -> list of anchor_alleles to check.
    mag_id : str
        The MAG ID we are processing, used to construct the profile file path.

    Returns:
    --------
    tuple : (sample_id, pd.Series)
        A tuple containing the sample_id and a pandas Series where the
        index is a MultiIndex (contig, position, anchor_allele) and
        the values are the frequencies. Returns (sample_id, None) on failure.
    """
    sample_id, sample_profile_dir = sample_data
    logger.debug(f"Starting processing for sample: {sample_id}")

    # Construct the profile path
    base_path = Path(sample_profile_dir) / f"{sample_id}_{mag_id}_profiled.tsv"
    # Use with_name to append .gz safely
    profile_path_gz = base_path.with_name(base_path.name + ".gz")

    if profile_path_gz.exists():
        profile_path = profile_path_gz
    elif base_path.exists():
        profile_path = base_path
    else:
        logger.warning(
            f"Profile file not found for {sample_id} and {mag_id}. "
            f"Checked: {profile_path_gz} and {base_path}"
        )
        return sample_id, None

    # Define dtypes for efficient profile file reading
    profile_dtypes = {
        "contig": str,
        "position": int,
        "total_coverage": float,
        **{nuc: "int32" for nuc in NUCLEOTIDES},
    }

    # Load the profile file for this sample
    profile_df = pd.read_csv(
        profile_path,
        sep="\t",
        dtype=profile_dtypes,
        usecols=list(profile_dtypes.keys()),
        compression="infer",  # Explicitly set to clarify intent
    )

    # Filter out sites with coverage below minimum threshold and set the index for fast lookup
    profile_df = profile_df[profile_df["total_coverage"] > min_cov_per_site - 1]
    profile_df = profile_df.set_index(["contig", "position"])
    if not profile_df.index.is_unique:
        raise ValueError(
            f"Profile file for {sample_id} and {mag_id} has duplicate (contig, position) entries."
        )
    # This dictionary will store the final frequency for each site
    results = {}

    # Iterate through all the sites we need to track
    for (contig, position), anchor_alleles in anchor_dict.items():
        try:
            # 1. Attempt to find the site
            site_row = profile_df.loc[(contig, position)]
            total_coverage = site_row["total_coverage"]

        except KeyError:
            # 2. Site is missing (e.g., 0 coverage). This is NOT an error.
            # Set all alleles for this site to NaN.
            logger.debug(
                f"Site {(contig, position)} not found in sample {sample_id}. Setting to NaN."
            )
            for allele in anchor_alleles:
                results[(contig, position, allele)] = np.nan
            continue  # Move to the next site in the loop

        # 3. Site was found. Now, process alleles.
        # This is intentionally NOT in the try-except block.
        # If an 'anchor_allele' is invalid (e.g., 'N'), we WANT this
        # to raise a KeyError and fail the worker, as per user's request.
        for allele in anchor_alleles:
            # This will raise KeyError if 'allele' is not in NUCLEOTIDES
            allele_count = site_row[allele]
            frequency = allele_count / total_coverage
            results[(contig, position, allele)] = frequency

    # Convert the results dict to a Series for easier combination later
    if not results:
        logger.warning(f"No anchor sites found in profile for {sample_id}")
        return sample_id, None

    result_series = pd.Series(results)
    result_series.name = sample_id
    logger.debug(f"Finished processing sample: {sample_id}")

    return sample_id, result_series


def process_samples_parallel(metadata_df, anchor_dict, mag_id, cpus, min_cov_per_site):
    """Handle multiprocessing logic for sample processing."""
    sample_data_list = [
        (row.sample_id, row.sample_profile_dir) for row in metadata_df.itertuples()
    ]

    # Use functools.partial to "freeze" the anchor_dict and mag_id arguments
    worker_func = partial(
        process_single_sample,
        anchor_dict=anchor_dict,
        mag_id=mag_id,
        min_cov_per_site=min_cov_per_site,
    )

    logger.info(f"Processing {len(sample_data_list)} samples across {cpus} CPUs...")

    # This will store the final Series from each successful worker
    all_results_dict = {}

    # Run Multiprocessing
    with multiprocessing.Pool(processes=cpus) as pool:
        # Use imap_unordered for progress tracking and efficiency
        for sample_id, result_series in tqdm(
            pool.imap_unordered(worker_func, sample_data_list),
            total=len(sample_data_list),
            desc="Processing Samples",
        ):
            if result_series is not None:
                all_results_dict[sample_id] = result_series

    return all_results_dict


def assemble_and_save_results(
    final_anchor_info_df,
    all_results_dict,
    output_dir,
    mag_id,
    anchor_column,
    metadata_df,
):
    """Combine results and save final output."""
    if not all_results_dict:
        logger.error("No samples were processed successfully. Exiting.")
        return

    logger.info("Sample processing complete. Assembling final table...")

    # Create the wide-form DataFrame
    freq_df = pd.DataFrame(all_results_dict)

    # Ensure index names are set for merging
    freq_df.index.names = ["contig", "position", "anchor_allele"]

    # Merge with the original anchor data to get gene_id, q_value, etc.
    final_df = final_anchor_info_df.merge(
        freq_df, left_index=True, right_index=True, how="left"
    )

    # Reset index to make contig/position/allele regular columns for R
    final_df = final_df.reset_index()

    # Add mag_id and anchor_criteria columns
    final_df["mag_id"] = mag_id
    final_df["anchor_criteria"] = anchor_column

    # Re-order columns for clarity
    id_cols = [
        "mag_id",
        "contig",
        "position",
        "gene_id",
        "anchor_allele",
        "anchor_criteria",
    ]

    # Add p-value columns if present
    if "min_p_value" in final_df.columns:
        id_cols.append("min_p_value")
    if "q_value" in final_df.columns:
        id_cols.append("q_value")

    # Get all other columns (sample_ids, etc.)
    other_cols = [col for col in final_df.columns if col not in id_cols]

    # Make sure gene_id is present, if not, don't add it
    if "gene_id" not in final_df.columns:
        id_cols.remove("gene_id")

    final_df = final_df[id_cols + other_cols]

    # Add QC column for number of samples with non-NaN values
    sample_cols = sorted(set(metadata_df["sample_id"]) & set(final_df.columns))
    final_df["n_samples_non_nan"] = final_df[sample_cols].notna().sum(axis=1)

    # Save the final wide-form table
    output_table_path = output_dir / f"{mag_id}_frequency_table.tsv"
    final_df.to_csv(output_table_path, sep="\t", index=False)

    logger.info(f"Success! Frequency table saved to: {output_table_path}")
    logger.info(
        f"Processed {len(final_anchor_info_df)} site-allele combinations "
        f"across {len(all_results_dict)} samples."
    )

    # Create and save the long-format version
    save_long_version(final_df, metadata_df, output_dir, mag_id)


def save_long_version(final_df, metadata_df, output_dir, mag_id):
    """Create a long-format version of the allele frequency tracking output.

    This function transforms the wide-format table (samples as columns) to long-format
    (one row per site-sample) and adds metadata columns for R/ggplot analysis.

    Parameters:
    -----------
    final_df : pd.DataFrame
        The final wide-format DataFrame with sample frequencies
    metadata_df : pd.DataFrame
        Metadata DataFrame with sample group/time/day information
    output_dir : Path
        Directory to save the long-format output
    mag_id : str
        MAG ID for naming the output file
    """
    id_cols = [
        "mag_id",
        "contig",
        "position",
        "gene_id",
        "anchor_allele",
        "anchor_criteria",
    ]

    # Add p-value columns if present
    if "min_p_value" in final_df.columns:
        id_cols.append("min_p_value")
    if "q_value" in final_df.columns:
        id_cols.append("q_value")

    id_cols = [c for c in id_cols if c in final_df.columns]
    sample_cols = sorted(set(metadata_df["sample_id"]) & set(final_df.columns))
    if not sample_cols:
        logger.warning("No sample columns found to melt; skipping long-format save.")
        return

    # value_vars = [c for c in final_df.columns if c not in id_cols]

    long_df = final_df.melt(
        id_vars=id_cols,
        value_vars=sample_cols,
        var_name="sample_id",
        value_name="frequency",
    )
    # attach group/time/subjectID and optional day/replicate
    merge_cols = ["sample_id", "group", "time", "subjectID"]
    if "day" in metadata_df.columns:
        merge_cols.append("day")
    if "replicate" in metadata_df.columns:
        merge_cols.append("replicate")
    long_df = long_df.merge(metadata_df[merge_cols], on="sample_id", how="left")
    long_path = output_dir / f"{mag_id}_frequency_table.long.tsv"
    long_df.to_csv(long_path, sep="\t", index=False)
    logger.info(f"Also wrote long format to: {long_path}")


def main():
    """Main function for the allele frequency tracking CLI."""
    parser = argparse.ArgumentParser(
        description=(
            "Track frequencies of pre-defined 'anchor alleles' for one MAG "
            "across all samples."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--mag-id",
        help="The MAG ID to process (e.g., 'SLG1007_DASTool_bins_61').",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--anchor-file",
        help=(
            "Path to the anchor allele file. This is the output from "
            "terminal_nucleotide_analysis.py (e.g., '..._terminal_nucleotides.tsv')."
        ),
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--metadata",
        help=(
            "Path to the enhanced metadata TSV file. Required columns: "
            "'sample_id', 'sample_profile_dir', 'group', 'time', 'subjectID'. "
            "Optional columns: 'day', 'replicate'."
        ),
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--output-dir",
        help="Directory to save the output files.",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--anchor-column",
        help="Column in the anchor-file to use as the anchor allele.",
        type=str,
        default="terminal_nucleotide_mean_freq",
        choices=[
            "terminal_nucleotide_mean_freq",
            "terminal_nucleotide_majority_vote",
            "terminal_nucleotide_freq_change",
        ],
    )
    parser.add_argument(
        "--cpus",
        help="Number of CPU cores to use for parallel processing.",
        type=int,
        default=multiprocessing.cpu_count(),
    )
    parser.add_argument(
        "--min-cov-per-site",
        help="Minimum total coverage per site to include in analysis. Sites with coverage below this threshold will be excluded.",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--log-level",
        help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    args = parser.parse_args()

    # Setup logging with the specified log level
    setup_logging(level=getattr(logging, args.log_level.upper()))

    logger.info("Starting allele frequency tracking workflow.")

    # Validate inputs and create output directory
    validate_inputs(args)

    # Load and validate metadata
    metadata_df = load_and_validate_metadata(args.metadata)

    # Load and prepare anchor sites
    final_anchor_info_df, anchor_dict = load_and_prepare_anchor_sites(
        args.anchor_file, args.mag_id, args.anchor_column
    )

    # Process samples in parallel
    all_results_dict = process_samples_parallel(
        metadata_df, anchor_dict, args.mag_id, args.cpus, args.min_cov_per_site
    )

    # Assemble and save results
    assemble_and_save_results(
        final_anchor_info_df,
        all_results_dict,
        args.output_dir,
        args.mag_id,
        args.anchor_column,
        metadata_df,
    )


if __name__ == "__main__":
    main()
