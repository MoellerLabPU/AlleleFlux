#!/usr/bin/env python3

"""
Positions-Based Quality Control Script for AlleleFlux

This script performs specialized quality control analysis on specific genomic positions
rather than full genome analysis. It's designed for targeted analysis of genomic regions
of interest such as significant SNPs, genes, or other genomic features.

Key Features:
- Position-specific breadth and coverage analysis
- Genome-wide baseline metrics for comparison
- Configurable denominator (positions vs genome size)
- Parallel processing for multiple samples
- Comprehensive logging and error handling
- Group-based aggregation of results

Key differences from main QC workflow (quality_control.py):
- ALWAYS requires a positions file (--positions_file parameter)
- Calculates BOTH positions-based AND genome-wide metrics for comparison
- DOES check genome-wide breadth/coverage thresholds (ensures overall sample quality)
- DOES aggregate results by experimental group
- NO timepoint validation (not meaningful for position subsets)
- NO subject/replicate counting (simplified workflow focused on position analysis)

Use Cases:
- Analyzing coverage at specific genomic positions of interest
- Quality control for variant calling or allele frequency analysis
- Comparing position-specific coverage patterns across samples
- Filtering samples based on genome-wide coverage baselines
"""

import argparse
import logging
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.qc_metrics import (
    aggregate_mag_stats,
    calculate_breadth_metrics,
    calculate_coverage_metrics,
    calculate_median_and_std_metrics,
    check_breadth_threshold,
    check_coverage_threshold,
    load_and_validate_profile,
)
from alleleflux.scripts.utilities.utilities import (
    build_contig_length_index,
    calculate_mag_sizes,
    load_mag_mapping,
    load_mag_metadata_file,
)

logger = logging.getLogger(__name__)


def load_positions_file(positions_file_path: Path) -> Dict[str, Set[Tuple[str, int]]]:
    """
    Load positions filter file and return MAG -> set of (contig, position) tuples.

    Parameters
    ----------
    positions_file_path : Path
        Path to TSV file with columns: MAG, contig, position

    Returns
    -------
    dict
        Mapping of MAG ID -> set of (contig, position) tuples

    Raises
    ------
    ValueError
        If duplicate positions are found for a MAG
    """
    pos_df = pd.read_csv(
        positions_file_path,
        sep="\t",
        dtype={"MAG": str, "contig": str, "position": "Int64"},
        usecols=["MAG", "contig", "position"],
    )

    # Check for duplicates
    dup_mask = pos_df.duplicated(subset=["MAG", "contig", "position"], keep=False)
    if dup_mask.any():
        dup_rows = pos_df[dup_mask].sort_values(["MAG", "contig", "position"])
        logger.error(f"Duplicate (MAG, contig, position) rows found:\n{dup_rows}")
        raise ValueError(
            f"Positions file contains {dup_mask.sum()} duplicate (MAG, contig, position) rows. "
            "Each position must appear only once per MAG."
        )

    # Build mapping MAG -> set of (contig, position) tuples for fast lookup
    positions_by_mag = {
        mag: set(
            zip(
                sub_df["contig"].astype(str),
                sub_df["position"].astype(int),
            )
        )
        for mag, sub_df in pos_df.groupby("MAG")
    }

    return positions_by_mag


def apply_positions_filter(
    df: pd.DataFrame,
    wanted_positions: Set[Tuple[str, int]],
    mag_id: str,
    sample_id: str,
) -> Tuple[pd.DataFrame, int]:
    """
    Filter dataframe to only rows in wanted_positions set.

    Parameters
    ----------
    df : pd.DataFrame
        Profile dataframe with 'contig' and 'position' columns
    wanted_positions : Set[Tuple[str, int]]
        Set of (contig, position) tuples to keep
    mag_id : str
        MAG identifier (for logging)
    sample_id : str
        Sample identifier (for logging)

    Returns
    -------
    tuple
        (filtered_df: pd.DataFrame, positions_universe_n: int)
        - filtered_df: DataFrame containing only wanted positions
        - positions_universe_n: Total number of positions in wanted set

    Raises
    ------
    ValueError
        If df lacks 'contig' or 'position' columns

    Notes
    -----
    - Converts contig to string and position to Int64 for consistent matching
    - Uses set membership for fast filtering
    - Logs info message if no matching positions found
    """
    if "contig" not in df.columns or "position" not in df.columns:
        raise ValueError(
            f"DataFrame must have 'contig' and 'position' columns for position filtering. "
            f"Available columns: {df.columns.tolist()}"
        )

    positions_universe_n = len(wanted_positions)

    # Ensure consistent data types for matching
    df["contig"] = df["contig"].astype(str)
    df["position"] = pd.to_numeric(df["position"]).astype("Int64")

    # Create boolean mask using set membership (much faster than merge)
    mask = df.apply(
        lambda row: (row["contig"], int(row["position"])) in wanted_positions, axis=1
    )
    filtered_df = df.loc[mask].copy()

    if filtered_df.empty:
        msg = (
            f"No matching (contig, position) rows for MAG {mag_id} in sample {sample_id} "
            f"after applying positions filter."
        )
        logger.info(msg)

    return filtered_df, positions_universe_n


def init_worker(
    metadata: dict,
    mag_sizes: dict,
    contig_to_mag_map: dict,
    positions_filter_map: dict,
    positions_denominator: str,
):
    """
    Initialize worker process with shared data.

    Parameters
    ----------
    metadata : dict
        Sample metadata dictionary
    mag_sizes : dict
        MAG ID -> size in bp
    contig_to_mag_map : dict
        Contig ID -> MAG ID mapping
    positions_filter_map : dict
        MAG ID -> set of (contig, position) tuples
    positions_denominator : str
        Either 'positions' or 'genome' for denominator choice
    """
    global metadata_dict
    global mag_size_dict
    global contig_to_mag
    global positions_filter
    global positions_den_mode

    metadata_dict = metadata
    mag_size_dict = mag_sizes
    contig_to_mag = contig_to_mag_map
    positions_filter = positions_filter_map
    positions_den_mode = positions_denominator


def process_positions_mag_sample(args):
    """
    Process a single sample for positions-based QC.

    This function calculates QC metrics for a specific set of positions while also
    capturing genome-wide metrics for comparison. It always uses a positions filter
    and calculates both filtered and genome-wide breadth/coverage.

    Parameters
    ----------
    args : tuple
        (sample_id, profile_path, mag_id, breadth_threshold, coverage_threshold)

    Returns
    -------
    dict
        Dictionary containing:
        - sample_id, MAG_ID, file_path, group, subjectID, replicate, time (if available)
        - genome_size: Size of MAG in bp
        - positions_considered: Number of positions in filter
        - breadth: Breadth for filtered positions (denominator depends on mode)
        - breadth_genome: Breadth across entire genome (always genome denominator)
        - breadth_threshold_passed: Based on genome-wide breadth
        - breadth_fail_reason: Explanation if failed
        - average_coverage: Average coverage for filtered positions
        - average_coverage_genome: Average coverage across entire genome
        - coverage_threshold_passed: Based on genome-wide coverage
        - coverage_fail_reason: Explanation if failed
        - median_coverage, median_coverage_including_zeros
        - coverage_std, coverage_std_including_zeros
    """
    sample_id, profile_path, mag_id, breadth_threshold, coverage_threshold = args

    # Build result dictionary
    # Columns are organized by category for better readability
    result = {
        # Sample metadata
        "sample_id": sample_id,
        "MAG_ID": mag_id,
        "file_path": profile_path,
        "group": metadata_dict[sample_id]["group"],
        "subjectID": metadata_dict[sample_id]["subjectID"],
        "replicate": metadata_dict[sample_id]["replicate"],
        "genome_size": mag_size_dict.get(mag_id, 0),
        "positions_considered": None,
        # Breadth metrics (grouped together)
        "breadth": None,
        "breadth_genome": None,
        "breadth_threshold_passed": False,
        "breadth_fail_reason": "",
        # Coverage metrics (grouped together)
        "average_coverage": None,
        "average_coverage_genome": None,
        "coverage_threshold_passed": False,
        "coverage_fail_reason": "",
        # Additional coverage statistics
        "median_coverage": None,
        "median_coverage_including_zeros": None,
        "coverage_std": None,
        "coverage_std_including_zeros": None,
    }

    # Add time if available
    if "time" in metadata_dict[sample_id]:
        result["time"] = metadata_dict[sample_id]["time"]

    # Check MAG size
    mag_size = mag_size_dict.get(mag_id)
    if mag_size is None or mag_size <= 0:
        msg = f"MAG size for {mag_id} not found or invalid."
        logger.error(f"{msg} sample={sample_id}, profile={profile_path}")
        result["breadth_fail_reason"] = msg
        return result

    # Get positions for this MAG
    wanted_positions = positions_filter.get(mag_id)
    if wanted_positions is None or len(wanted_positions) == 0:
        msg = f"No positions found in filter for MAG {mag_id}"
        logger.error(msg)
        result["breadth_fail_reason"] = msg
        return result

    # Calculate positions_considered count for this MAG (same for all samples of this MAG)
    positions_considered_count = len(wanted_positions)

    # Load and validate profile
    try:
        df = load_and_validate_profile(profile_path, mag_id, sample_id, contig_to_mag)
    except ValueError as e:
        result["breadth_fail_reason"] = str(e)
        return result

    # Check for zero coverage values in the profile (should not exist)
    zero_coverage_count = (df["total_coverage"] == 0).sum()
    if zero_coverage_count > 0:
        logger.warning(
            f"Profile contains {zero_coverage_count} positions with zero coverage for MAG {mag_id}. "
            f"These should typically be absent from the profile. file={profile_path} sample={sample_id}"
        )
    positions_with_coverage_genome = df[df["total_coverage"] >= 1].shape[0]
    total_coverage_sum_genome = df["total_coverage"].sum()

    # Apply positions filter
    filtered_df, positions_universe_n = apply_positions_filter(
        df, wanted_positions, mag_id, sample_id
    )

    if filtered_df.empty:
        result["positions_considered"] = positions_considered_count
        msg = f"No matching positions found for MAG {mag_id} in sample {sample_id}"
        logger.info(msg)
        result["breadth_fail_reason"] = msg
        return result

    result["positions_considered"] = positions_considered_count

    # Determine denominator based on mode
    if positions_den_mode == "positions":
        denom = positions_universe_n
    else:  # genome mode
        denom = mag_size

    # Calculate breadth metrics (filtered + genome-wide)
    breadth_metrics = calculate_breadth_metrics(
        filtered_df, denom, mag_size, True, positions_with_coverage_genome
    )
    result.update(breadth_metrics)

    # Calculate coverage metrics (filtered + genome-wide)
    coverage_metrics = calculate_coverage_metrics(
        filtered_df, denom, mag_size, True, total_coverage_sum_genome
    )
    result.update(coverage_metrics)

    # Calculate median and std metrics (for filtered positions)
    avg_cov = result["average_coverage"]
    median_std_metrics = calculate_median_and_std_metrics(filtered_df, denom, avg_cov)
    result.update(median_std_metrics)

    # Check genome-wide thresholds
    # When using position filter: ALWAYS check genome-wide metrics (breadth_genome, avg_cov_genome)
    # This ensures samples have adequate genome-wide coverage even for position subsets
    breadth_passed, breadth_fail_reason = check_breadth_threshold(
        result["breadth"],
        result["breadth_genome"],
        breadth_threshold,
        using_filter=True,
        sample_id=sample_id,
        mag_id=mag_id,
    )
    result["breadth_threshold_passed"] = breadth_passed
    result["breadth_fail_reason"] = breadth_fail_reason

    # Check coverage threshold (always uses genome-wide metric when filtering)
    coverage_passed, coverage_fail_reason = check_coverage_threshold(
        result["average_coverage"],
        result["average_coverage_genome"],
        coverage_threshold,
        breadth_passed,
        using_filter=True,
        sample_id=sample_id,
        mag_id=mag_id,
    )
    result["coverage_threshold_passed"] = coverage_passed
    result["coverage_fail_reason"] = coverage_fail_reason

    return result


def process_mag(args):
    """
    Process all samples for a single MAG.

    Parameters
    ----------
    args : tuple
        Contains all parameters needed for processing

    Returns
    -------
    tuple
        (group_summary_df, group_time_summary_df, overall_summary_df)
    """
    (
        mag_id,
        metadata_file,
        mag_size_dict,
        contig_to_mag_map,
        positions_filter_map,
        output_dir,
        breadth_threshold,
        coverage_threshold,
        cpus,
        positions_denominator,
        data_type,
    ) = args

    logger.info(f"Processing MAG: {mag_id}")

    # Load metadata - use "single" data_type since we don't do timepoint validation
    metadata_dict, sample_files_with_mag_id = load_mag_metadata_file(
        metadata_file, mag_id, breadth_threshold, coverage_threshold, data_type
    )

    if not sample_files_with_mag_id:
        raise ValueError(
            f"No samples found in metadata file {metadata_file} for MAG {mag_id}"
        )

    num_procs = min(cpus, len(sample_files_with_mag_id))
    logger.info(
        f"Processing {len(sample_files_with_mag_id)} samples for {mag_id} "
        f"with {num_procs} processes"
    )

    # Process samples in parallel
    with Pool(
        processes=num_procs,
        initializer=init_worker,
        initargs=(
            metadata_dict,
            mag_size_dict,
            contig_to_mag_map,
            positions_filter_map,
            positions_denominator,
        ),
    ) as pool:
        results_list = list(
            pool.imap_unordered(process_positions_mag_sample, sample_files_with_mag_id)
        )

    # Build results DataFrame
    df_results = pd.DataFrame(results_list)

    # Save sample-level results
    out_file = Path(output_dir) / f"{mag_id}_positions_QC.tsv"
    df_results.to_csv(out_file, sep="\t", index=False)
    logger.info(f"Saved positions QC report for {mag_id} to {out_file}")

    # Aggregate summaries (by group and overall)
    group_summary = aggregate_mag_stats(
        df_results, mag_id, group_cols=["group"], overall=False
    )

    group_time_summary = None
    if "time" in df_results.columns:
        group_time_summary = aggregate_mag_stats(
            df_results, mag_id, group_cols=["group", "time"], overall=False
        )

    overall_summary = aggregate_mag_stats(
        df_results, mag_id, group_cols=[], overall=True
    )

    return group_summary, group_time_summary, overall_summary


def main():
    """Main entry point for positions-based QC script."""
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Perform quality control on specific genomic positions for MAGs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--metadata_dir",
        required=True,
        type=Path,
        help="Directory containing *_metadata.tsv files (one per MAG)",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="FASTA file with contigs (for MAG sizes)",
    )
    parser.add_argument(
        "--mag_mapping_file",
        required=True,
        type=Path,
        help="Tab-separated file mapping contig names to MAG IDs (columns: contig_id, mag_id)",
    )
    parser.add_argument(
        "--positions_file",
        required=True,
        type=Path,
        help="TSV file with columns: MAG, contig, position. Specifies which positions to analyze.",
    )
    parser.add_argument(
        "--positions_denominator",
        choices=["positions", "genome"],
        default="positions",
        help=(
            "Denominator for breadth/coverage calculations: 'positions' uses the number of "
            "specified positions, 'genome' uses the full genome size"
        ),
    )
    parser.add_argument(
        "--data_type",
        help="Is the data from a single timepoint or from a time series (longitudinal)",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )
    parser.add_argument(
        "--breadth_threshold",
        type=float,
        default=0.1,
        help="Minimum genome-wide breadth coverage required to pass (0.0-1.0)",
    )
    parser.add_argument(
        "--coverage_threshold",
        type=float,
        default=1.0,
        help="Minimum genome-wide average coverage depth required to pass",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=cpu_count(),
        help="Number of processors to use",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="Directory to write output files",
    )

    args = parser.parse_args()

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Calculate MAG sizes and build contig indices
    logger.info("Loading MAG information...")
    mag_size_dict = calculate_mag_sizes(args.fasta, args.mag_mapping_file)
    contig_to_mag_global = load_mag_mapping(args.mag_mapping_file)
    # Note: contig_length_dict not needed for positions-based QC (no length-weighted coverage calculation)

    # Load positions filter (required for this script)
    logger.info(f"Loading positions from {args.positions_file}...")
    positions_by_mag = load_positions_file(args.positions_file)
    logger.info(
        f"Loaded positions for {len(positions_by_mag)} MAG(s): "
        f"{sum(len(v) for v in positions_by_mag.values())} total positions"
    )

    # Find metadata files
    metadata_files = list(args.metadata_dir.glob("*_metadata.tsv"))
    if not metadata_files:
        raise ValueError(f"No *_metadata.tsv files found in {args.metadata_dir}")

    # Filter to only MAGs that have positions in the filter
    mags_in_positions = set(positions_by_mag.keys())
    metadata_files = [
        mf
        for mf in metadata_files
        if mf.stem.split("_metadata")[0] in mags_in_positions
    ]

    if not metadata_files:
        raise ValueError(
            f"No MAGs from {args.metadata_dir} found in positions file. "
            "Check that MAG IDs in positions file match metadata filenames."
        )

    logger.info(f"Processing {len(metadata_files)} MAG(s) with position filters")

    # Build tasks for each MAG
    tasks = [
        (
            metadata_file.stem.split("_metadata")[0],  # mag_id from filename
            metadata_file,
            mag_size_dict,
            contig_to_mag_global,
            positions_by_mag,
            args.output_dir,
            args.breadth_threshold,
            args.coverage_threshold,
            args.cpus,
            args.positions_denominator,
            args.data_type,
        )
        for metadata_file in metadata_files
    ]

    # Process each MAG and collect summaries
    all_group = []
    all_group_time = []
    all_overall = []

    for task in tasks:
        group_df, group_time_df, overall_df = process_mag(task)
        if group_df is not None:
            all_group.append(group_df)
        if group_time_df is not None:
            all_group_time.append(group_time_df)
        if overall_df is not None:
            all_overall.append(overall_df)

    # Write combined summary files
    summaries = [
        (all_group, "ALL_MAGs_positions_QC_group_summary.tsv"),
        (all_group_time, "ALL_MAGs_positions_QC_group_time_summary.tsv"),
        (all_overall, "ALL_MAGs_positions_QC_overall_summary.tsv"),
    ]

    for dfs, filename in summaries:
        if dfs:
            combined_df = pd.concat(dfs, ignore_index=True)
            output_path = args.output_dir / filename
            combined_df.to_csv(output_path, sep="\t", index=False)
            logger.info(f"Saved combined summary to {output_path}")

    logger.info("Positions QC analysis completed successfully.")


if __name__ == "__main__":
    main()
