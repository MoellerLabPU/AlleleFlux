#!/usr/bin/env python3
"""
Terminal nucleotide analysis for specified group and timepoint samples.

This script processes genomic data to identify terminal nucleotides in specified samples
by analyzing nucleotide frequency profiles for significant sites.

Enhanced with major allele majority voting algorithm alongside the existing mean frequency approach.
"""

import argparse
import logging
import math
import multiprocessing
from collections import Counter
from functools import partial
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# Global constant for nucleotide symbols
NUCLEOTIDES = ["A", "C", "G", "T"]


def find_profile_files(profile_dir, sample_ids, mag_ids=None):
    """
    Find profile files in subdirectories by sample ID and MAG ID.

    Parameters:
    -----------
    profile_dir : str
        Root directory containing profile subdirectories by sample ID
    sample_ids : list
        List of sample IDs to look for
    mag_ids : list, optional
        List of MAG IDs to filter profile files by. If None, uses wildcard pattern.

    Returns:
    --------
    dict : Dictionary mapping sample_id to list of profile file paths
    """
    profile_files = {}
    profile_dir_path = Path(profile_dir)

    if not profile_dir_path.exists():
        raise FileNotFoundError(f"Profile directory not found: {profile_dir}")

    logger.info(
        f"Searching for profile files for {len(sample_ids)} sample(s) and {len(mag_ids)} MAG(s)"
    )

    # Check if profile_dir contains subdirectories (new structure) or files directly
    subdirs = [d for d in profile_dir_path.iterdir() if d.is_dir()]

    if subdirs:
        # New structure: subdirectories by sample_id
        logger.info(f"Found {len(subdirs)} subdirectories in profile directory")
        for subdir in subdirs:
            sample_id = subdir.name
            if sample_id in sample_ids:
                for mag_id in mag_ids:
                    base = f"{sample_id}_{mag_id}_profiled.tsv"
                    candidates = list(subdir.glob(base)) + list(
                        subdir.glob(base + ".gz")
                    )
                    if candidates:
                        profile_files.setdefault(sample_id, []).extend(
                            map(str, candidates)
                        )
    # Log summary
    total_files_found = sum(len(files) for files in profile_files.values())
    samples_with_files = len(profile_files)

    logger.info(
        f"Profile file search complete: found {total_files_found} files for {samples_with_files}/{len(sample_ids)} samples"
    )

    return profile_files


def get_terminal_samples(metadata_file, group, timepoint):
    """
    Get sample IDs that match the specified group and timepoint.

    Parameters:
    -----------
    metadata_file : str
        Path to metadata file
    group : str
        Group name to filter for
    timepoint : str
        Timepoint to filter for

    Returns:
    --------
    list : List of sample IDs matching the criteria
    """
    logger.info(
        f"Loading metadata to identify samples for group '{group}' and timepoint '{timepoint}'"
    )
    metadata = pd.read_csv(metadata_file, sep="\t")

    if (
        "group" not in metadata.columns
        or "time" not in metadata.columns
        or "sample_id" not in metadata.columns
    ):
        raise ValueError(
            "Metadata must contain 'group', 'time', and 'sample_id' columns"
        )

    # Filter for matching group and timepoint
    matching_samples = metadata[
        (metadata["group"] == group) & (metadata["time"] == timepoint)
    ]["sample_id"].tolist()

    if not matching_samples:
        raise ValueError(
            f"No samples found for group '{group}' and timepoint '{timepoint}' in metadata"
        )

    logger.info(
        f"Found {len(matching_samples)} matching samples: {matching_samples[:5]}{'...' if len(matching_samples) > 5 else ''}"
    )

    return matching_samples


def perform_majority_voting(major_alleles):
    """
    Perform majority voting to determine terminal nucleotide from major alleles.

    When a sample has tied major alleles, all tied alleles contribute fractionally
    to the vote count (each gets 1/n_tied votes) to avoid biasing toward any single nucleotide.

    Parameters:
    -----------
    major_alleles : list
        List of major allele strings from all samples, or None values for missing data.
        Each string can be a single nucleotide or comma-separated tied nucleotides.

    Returns:
    --------
    str : The terminal nucleotide(s) determined by majority voting (comma-separated if tied)
    """
    # Count votes with fractional contribution for ties
    allele_counts = Counter()

    for major_allele in major_alleles:
        if major_allele is None:
            continue
        # major_allele can be single or comma-separated nucleotides
        votes = major_allele.split(",")
        # Each tied nucleotide gets equal fractional vote
        vote_weight = 1.0 / len(votes)

        for nucleotide in votes:
            allele_counts[nucleotide] += vote_weight

    # Handle "no votes at all"
    if not allele_counts:
        return None

    # Find maximum vote count
    max_count = max(allele_counts.values())

    # Get all nucleotides with maximum count and sort for deterministic output
    tied_nucleotides = sorted(
        [
            nuc
            for nuc, count in allele_counts.items()
            if math.isclose(count, max_count, rel_tol=1e-9)
        ]
    )

    if len(tied_nucleotides) > 1:
        logger.debug(
            f"Tie detected in majority voting: {tied_nucleotides} with count {max_count:.2f}"
        )

    # Return all tied nucleotides as comma-separated string
    return ",".join(tied_nucleotides)


def worker_wrapper(mag_data, profile_files, sample_ids, p_value_column, output_dir):
    """
    Worker function wrapper to unpack tuple arguments for parallel processing.

    This must be a top-level function to be picklable on macOS/Windows (spawn start method).

    Parameters:
    -----------
    mag_data : tuple
        Tuple of (mag_id, mag_sites_dataframe)
    profile_files : dict
        Dictionary mapping sample_id to profile file paths
    sample_ids : list
        List of sample IDs to process
    p_value_column : str
        Column name for p-values
    output_dir : str
        Output directory

    Returns:
    --------
    dict : Results summary from process_single_mag
    """
    mag_id, mag_sites = mag_data
    return process_single_mag(
        mag_id=mag_id,
        mag_sites=mag_sites,
        profile_files=profile_files,
        sample_ids=sample_ids,
        p_value_column=p_value_column,
        output_dir=output_dir,
    )


def process_single_mag(
    mag_id,
    mag_sites,
    profile_files,
    sample_ids,
    p_value_column,
    output_dir,
):
    """
    Process a single MAG's significant sites.

    Parameters:
    -----------
    mag_id : str
        MAG identifier
    mag_sites : pd.DataFrame
        Sites for this MAG
    profile_files : dict
        Dictionary mapping sample_id to profile file paths
    sample_ids : list
        List of sample IDs to process
    p_value_column : str
        Column name for p-values
    output_dir : str
        Output directory

    Returns:
    --------
    dict : Results summary
    """
    logger.info(f"Processing MAG: {mag_id} with {len(mag_sites)} sites")

    # Define dtypes for profile files
    profile_dtypes = {
        "contig": str,
        "position": int,
        "total_coverage": float,
        **{nuc: "int32" for nuc in NUCLEOTIDES},
    }

    # Cache profile data per sample to avoid repeated file I/O
    sample_profile_cache = {}
    for sample_id in sample_ids:
        if sample_id not in profile_files:
            logger.warning(f"No profile files found for sample {sample_id}")
            continue

        # Filter profile files for this specific MAG ID
        # Profile files are named: {sample_id}_{mag_id}_profiled.tsv or .tsv.gz
        mag_specific_files = [
            f for f in profile_files[sample_id] if f"_{mag_id}_profiled.tsv" in f
        ]

        if len(mag_specific_files) == 0:
            logger.warning(
                f"No profile file found for sample {sample_id} and MAG {mag_id}"
            )
            continue
        elif len(mag_specific_files) > 1:
            raise ValueError(
                f"Multiple profile files found for sample {sample_id} and MAG {mag_id}: {mag_specific_files}. Using first match."
            )

        profile_file = mag_specific_files[0]
        profile_data = pd.read_csv(
            profile_file,
            sep="\t",
            dtype=profile_dtypes,
            usecols=profile_dtypes.keys(),
        )
        # Filter for positions with total_coverage > 0
        profile_data = profile_data[profile_data["total_coverage"] > 0]
        # Index by (contig, position) for fast lookup
        profile_data = profile_data.set_index(["contig", "position"])
        sample_profile_cache[sample_id] = profile_data

    # Process each significant site for this MAG
    terminal_nucleotides_mean_freq = []
    majority_vote_terminal = []
    intermediate_data = []

    for idx, site in mag_sites.iterrows():
        contig = site["contig"]
        position = site["position"]

        # Collect nucleotide frequencies for this site across all samples
        site_frequencies = {}
        sample_major_alleles = []
        samples_with_site = []  # Track which samples have this site present

        for sample_id in sample_ids:
            if sample_id not in sample_profile_cache:
                # Sample data not available - exclude from mean calculation
                sample_major_alleles.append(None)
                continue

            # Fast lookup using multi-index
            profile_data = sample_profile_cache[sample_id]
            try:
                site_row = profile_data.loc[(contig, position)]

                total_coverage = site_row["total_coverage"]

                # Calculate frequencies using vectorized operation
                sample_nucleotide_frequencies = {
                    nuc: site_row[nuc] / total_coverage for nuc in NUCLEOTIDES
                }

                # Store in site_frequencies dict
                for nuc, freq in sample_nucleotide_frequencies.items():
                    site_frequencies[f"{nuc}_{sample_id}"] = freq

                # Identify major allele(s) - include all nucleotides with maximum frequency
                max_freq = max(sample_nucleotide_frequencies.values())
                major_alleles = sorted(
                    [
                        nuc
                        for nuc, freq in sample_nucleotide_frequencies.items()
                        if math.isclose(freq, max_freq, rel_tol=1e-9)
                    ]
                )
                # Store as comma-separated string (for both output and voting)
                # When there are ties, all tied alleles are passed to voting with equal weight
                major_allele = ",".join(major_alleles)
                sample_major_alleles.append(major_allele)
                samples_with_site.append(sample_id)  # Site is present in this sample
            except KeyError:
                # Site not found in this profile - exclude from mean calculation
                logger.debug(
                    f"Site {contig}:{position} not found in profile for sample {sample_id}. Excluding from mean calculation."
                )
                sample_major_alleles.append(None)

        # Calculate mean frequencies using only samples where the site is present
        if samples_with_site:
            mean_frequencies = {
                nuc: sum(
                    site_frequencies.get(f"{nuc}_{sample_id}", 0)
                    for sample_id in samples_with_site
                )
                / len(samples_with_site)
                for nuc in NUCLEOTIDES
            }
        else:
            # No samples have this site - set all frequencies to NaN to distinguish from observed zeros
            mean_frequencies = {nuc: math.nan for nuc in NUCLEOTIDES}

        # Store intermediate data - one row per site using dict comprehension
        site_row = {
            "mag_id": mag_id,
            "contig": contig,
            "position": position,
            "gene_id": site["gene_id"],
            p_value_column: site[p_value_column],
            # Add all nucleotide-sample frequency columns
            **{
                f"{nuc}_{sample_id}": site_frequencies.get(
                    f"{nuc}_{sample_id}", math.nan
                )
                for nuc in NUCLEOTIDES
                for sample_id in sample_ids
            },
            # Add mean frequency columns
            **{f"{nuc}_mean_frequency": mean_frequencies[nuc] for nuc in NUCLEOTIDES},
            "n_samples_used_for_mean": len(samples_with_site),
            # Add major allele columns
            **{
                f"major_allele_{sample_id}": (
                    sample_major_alleles[i]
                    if i < len(sample_major_alleles)
                    and sample_major_alleles[i] is not None
                    else None
                )
                for i, sample_id in enumerate(sample_ids)
            },
        }
        # Method 1: Terminal nucleotide based on highest mean frequency
        if not samples_with_site:
            mean_freq_terminal_nuc = None
        else:
            max_mean_freq = max(mean_frequencies.values())
            terminal_nucleotides_tied = sorted(
                [
                    nuc
                    for nuc, freq in mean_frequencies.items()
                    if math.isclose(freq, max_mean_freq, rel_tol=1e-9)
                ]
            )
            mean_freq_terminal_nuc = (
                ",".join(terminal_nucleotides_tied)
                if len(terminal_nucleotides_tied) > 1
                else terminal_nucleotides_tied[0]
            )

        # Method 2: Majority voting approach
        majority_vote_terminal_nuc = perform_majority_voting(sample_major_alleles)

        site_row["terminal_nucleotide_mean_freq"] = mean_freq_terminal_nuc
        site_row["terminal_nucleotide_majority_vote"] = majority_vote_terminal_nuc

        intermediate_data.append(site_row)
        terminal_nucleotides_mean_freq.append(mean_freq_terminal_nuc)
        majority_vote_terminal.append(majority_vote_terminal_nuc)

    # Create output DataFrames
    processed_sites = mag_sites.copy()
    processed_sites["terminal_nucleotide_mean_freq"] = terminal_nucleotides_mean_freq
    processed_sites["terminal_nucleotide_majority_vote"] = majority_vote_terminal

    # Drop source_file column if it exists
    if "source_file" in processed_sites.columns:
        processed_sites = processed_sites.drop(columns=["source_file"])

    intermediate_df = pd.DataFrame(intermediate_data)

    # Create output directories
    mag_output_dir = Path(output_dir) / mag_id
    mag_output_dir.mkdir(parents=True, exist_ok=True)

    # Save results
    main_output_file = mag_output_dir / f"{mag_id}_terminal_nucleotides.tsv"
    processed_sites.to_csv(main_output_file, sep="\t", index=False)

    intermediate_output_file = (
        mag_output_dir / f"{mag_id}_nucleotide_frequencies.tsv.gz"
    )
    intermediate_df.to_csv(
        intermediate_output_file, sep="\t", index=False, compression="gzip"
    )

    # Generate summary statistics (per-nucleotide counts, splitting ties)
    def _split_count(series: pd.Series) -> dict:
        counter: dict[str, int] = {}
        for val in series.dropna():
            for nuc in str(val).split(","):
                if nuc == "" or nuc not in NUCLEOTIDES:
                    continue
                counter[nuc] = counter.get(nuc, 0) + 1
        return counter

    terminal_counts = _split_count(processed_sites["terminal_nucleotide_mean_freq"])
    majority_counts = _split_count(processed_sites["terminal_nucleotide_majority_vote"])

    logger.info(f"MAG {mag_id} complete.")
    logger.info(f"  Mean frequency method: {terminal_counts}")
    logger.info(f"  Majority voting method: {majority_counts}")

    return {
        "mag_id": mag_id,
        "sites_processed": len(processed_sites),
        "mean_freq_terminal_nucleotides": terminal_counts,
        "majority_vote_terminal_nucleotides": majority_counts,
        "main_output": str(main_output_file),
        "intermediate_output": (
            str(intermediate_output_file) if not intermediate_df.empty else None
        ),
        "samples_processed": len(sample_ids),
    }


def main():
    """Main function for terminal nucleotide analysis CLI."""
    # Parse arguments first to get log level
    parser = argparse.ArgumentParser(
        description="Identify terminal nucleotides in specified group and timepoint samples using dual methodology",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--significant_sites",
        help="Tab-separated file containing significant sites with columns 'mag_id', 'contig', 'position', 'gene_id', and either 'min_p_value' or 'min_q_value'",
        type=str,
        required=True,
        metavar="filepath",
    )

    parser.add_argument(
        "--profile_dir",
        help="Directory containing nucleotide frequency profiles for each sample (organized in subdirectories by sample_id)",
        type=str,
        required=True,
        metavar="directory",
    )

    parser.add_argument(
        "--group",
        help="Specific group name to identify terminal timepoint samples",
        type=str,
        required=True,
        metavar="string",
    )

    parser.add_argument(
        "--timepoint",
        help="Specific timepoint to identify terminal samples",
        type=str,
        required=True,
        metavar="string",
    )

    parser.add_argument(
        "--metadata",
        help="Path to metadata table (same format as mag_metadata.py) containing sample information",
        type=str,
        required=True,
        metavar="filepath",
    )

    parser.add_argument(
        "--output",
        help="Output directory for results",
        type=str,
        required=True,
        metavar="directory",
    )

    parser.add_argument(
        "--p_value_column",
        help="Column to use for p-value filtering",
        type=str,
        choices=["min_p_value", "q_value"],
        default="q_value",
        metavar="string",
    )

    parser.add_argument(
        "--p_value_threshold",
        type=float,
        default=0.05,
        help="Threshold for significance filtering using the selected p-value column",
        metavar="float",
    )

    parser.add_argument(
        "--test-type",
        type=str,
        choices=[
            "two_sample_unpaired_tTest",
            "two_sample_unpaired_MannWhitney",
            "two_sample_unpaired_tTest_abs",
            "two_sample_unpaired_MannWhitney_abs",
            "two_sample_paired_tTest",
            "two_sample_paired_Wilcoxon",
            "two_sample_paired_tTest_abs",
            "two_sample_paired_Wilcoxon_abs",
            "single_sample_tTest",
            "single_sample_Wilcoxon",
            "CMH",
            "LMM",
            "LMM_abs",
            # "lmm_across_time",
            # "cmh_across_time",
        ],
        default="two_sample_paired_tTest",
        help="A test type to filter the `significant_sites` dataframe.",
    )

    parser.add_argument(
        "--group_analyzed",
        help="Optional group name to filter significant sites by the 'group_analyzed' column if present",
        type=str,
        default=None,
        metavar="string",
    )

    parser.add_argument(
        "--cpus",
        help="Number of CPU cores to use for parallel processing (default: all available)",
        type=int,
        default=multiprocessing.cpu_count(),
        metavar="integer",
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level",
        metavar="string",
    )

    args = parser.parse_args()

    # Setup logging with the specified level
    setup_logging(level=getattr(logging, args.log_level))

    # Validate input files exist
    input_files = [args.significant_sites, args.metadata]
    for file_path in input_files:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

    if not Path(args.profile_dir).exists():
        raise FileNotFoundError(f"Profile directory not found: {args.profile_dir}")

    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("Starting enhanced terminal nucleotide analysis")

    # Load significant sites
    logger.info("Loading significant sites table")
    significant_sites = pd.read_csv(args.significant_sites, sep="\t")

    if significant_sites.empty:
        raise ValueError("Significant sites file is empty")

    # Validate required columns
    required_cols = ["mag_id", "contig", "position", "gene_id", args.p_value_column]
    missing_cols = set(required_cols) - set(significant_sites.columns)
    if missing_cols:
        raise ValueError(
            f"Missing required columns in significant sites: {missing_cols}"
        )

    # Apply significance threshold filtering
    logger.info(
        f"Applying significance threshold: {args.p_value_column} <= {args.p_value_threshold}"
    )
    initial_count = len(significant_sites)
    significant_sites = significant_sites[
        significant_sites[args.p_value_column] <= args.p_value_threshold
    ].copy()
    logger.info(
        f"After significance filtering ({args.p_value_column} <= {args.p_value_threshold}), "
        f"{len(significant_sites)}/{initial_count} sites remain"
    )

    if significant_sites.empty:
        raise ValueError(
            f"No sites remain after applying {args.p_value_column} <= {args.p_value_threshold} threshold"
        )

    # Filter by test_type if specified
    if args.test_type:
        if "test_type" in significant_sites.columns:
            initial_count = len(significant_sites)
            significant_sites = significant_sites[
                significant_sites["test_type"] == args.test_type
            ].copy()
            logger.info(
                f"Filtered by test_type '{args.test_type}': {len(significant_sites)}/{initial_count} sites remain"
            )
            if significant_sites.empty:
                raise ValueError(f"No sites found for test_type '{args.test_type}'")
        else:
            logger.warning(
                f"test_type column not found in significant_sites file. Ignoring --test_type filter."
            )

    # Filter by group_analyzed if specified
    if args.group_analyzed:
        if "group_analyzed" in significant_sites.columns:
            initial_count = len(significant_sites)
            significant_sites = significant_sites[
                significant_sites["group_analyzed"] == args.group_analyzed
            ].copy()
            logger.info(
                f"Filtered by group_analyzed '{args.group_analyzed}': {len(significant_sites)}/{initial_count} sites remain"
            )
            if significant_sites.empty:
                raise ValueError(
                    f"No sites found for group_analyzed '{args.group_analyzed}'"
                )
        else:
            logger.warning(
                f"group_analyzed column not found in significant_sites file. Ignoring --group_analyzed filter."
            )

    # Get matching samples from metadata
    terminal_samples = get_terminal_samples(args.metadata, args.group, args.timepoint)

    # Find profile files for specific MAG IDs from significant sites
    mag_groups = significant_sites.groupby("mag_id")
    mag_ids = list(mag_groups.groups.keys())
    profile_files = find_profile_files(args.profile_dir, terminal_samples, mag_ids)

    if not profile_files:
        raise ValueError(
            f"No profile files found for terminal samples in {args.profile_dir}"
        )

    # Optimize CPU usage: use minimum of specified CPUs and number of MAGs
    # to avoid wasting resources when there are fewer MAGs than CPUs
    effective_cpus = min(args.cpus, len(mag_ids))

    logger.info(
        f"Processing {len(mag_ids)} MAG(s) with significant sites using {effective_cpus} CPU(s)"
    )

    # Prepare MAG data for parallel processing
    # Create list of tuples: [(mag_id, mag_sites_dataframe), (mag_id, mag_sites_dataframe), ...]
    # Each tuple contains the MAG ID and its corresponding DataFrame of significant sites
    mag_data_list = [(mag_id, mag_groups.get_group(mag_id)) for mag_id in mag_ids]

    # Process MAGs in parallel using multiprocessing
    all_results = []
    with multiprocessing.Pool(processes=effective_cpus) as pool:
        # Use functools.partial to bind shared arguments to top-level worker function
        # This ensures the worker is picklable on macOS/Windows (spawn start method)
        func = partial(
            worker_wrapper,
            profile_files=profile_files,
            sample_ids=terminal_samples,
            p_value_column=args.p_value_column,
            output_dir=args.output,
        )

        # Execute parallel processing with streaming results using imap_unordered
        # Results are yielded as they complete rather than waiting for all to finish
        for result in tqdm(
            pool.imap_unordered(func, mag_data_list),
            total=len(mag_data_list),
            desc="Processing MAGs",
        ):
            # Filter out None results (failed MAGs) and collect successful ones
            if result is not None:
                all_results.append(result)

    # Create summary report
    if all_results:
        summary_data = []
        for result in all_results:
            summary_data.append(
                {
                    "mag_id": result["mag_id"],
                    "sites_processed": result["sites_processed"],
                    "samples_processed": result["samples_processed"],
                    "main_output": result["main_output"],
                    "intermediate_output": result["intermediate_output"],
                    **{
                        f"mean_freq_terminal_{nuc}_count": count
                        for nuc, count in result[
                            "mean_freq_terminal_nucleotides"
                        ].items()
                    },
                    **{
                        f"majority_vote_terminal_{nuc}_count": count
                        for nuc, count in result[
                            "majority_vote_terminal_nucleotides"
                        ].items()
                    },
                }
            )

        summary_df = pd.DataFrame(summary_data)
        summary_file = output_path / "terminal_nucleotide_analysis_summary.tsv"
        summary_df.to_csv(summary_file, sep="\t", index=False)

        logger.info(f"Analysis complete. Processed {len(all_results)} MAG(s)")
        logger.info(f"Summary report saved to: {summary_file}")

        # Log overall statistics
        total_sites = sum(r["sites_processed"] for r in all_results)
        total_mags = len(all_results)
        logger.info(f"Total: {total_sites} sites across {total_mags} MAG(s)")
        logger.info(
            "Both mean frequency and majority voting methods completed successfully"
        )

    else:
        logger.warning("No MAGs were successfully processed")


if __name__ == "__main__":
    main()
