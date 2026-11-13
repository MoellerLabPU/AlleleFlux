#!/usr/bin/env python3

"""
Genome-wide Quality Control Script for AlleleFlux

This script performs comprehensive quality control analysis on metagenome-assembled genomes (MAGs).
It calculates coverage statistics, validates sample quality, and performs threshold-based filtering.

Key Features:
- Genome-wide breadth and coverage analysis
- Length-weighted coverage calculations
- Parallel processing for multiple samples
- Comprehensive logging and error handling

The script processes coverage profile files for each MAG and determines whether samples pass
breadth and coverage thresholds. It supports both single timepoint and longitudinal study designs.

The tables it writes hold PER-SAMPLE facts only (one row per sample of the
metadata folder it was given, normally the samples of one timepoint).  Which
subjects are readable at both timepoints of a comparison, and how many
replicates each group then has, are questions about a comparison and its two
tested groups: ``eligibility_table.py`` answers them from these tables.

For position-specific QC analysis (analyzing coverage at specific genomic positions),
use the separate positions_qc.py script in alleleflux/scripts/accessory/.

Author: AlleleFlux Development Team
"""

import argparse
import logging
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.qc_metrics import (
    aggregate_mag_stats,
    calculate_breadth_metrics,
    calculate_coverage_metrics,
    calculate_length_weighted_coverage,
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

# The columns of a per-genome QC table ({mag}_QC.tsv), in the order they are
# written: who the sample is, then breadth, then coverage.  Longitudinal data
# add ``time`` at the end (qc_table_columns).  The list exists so that a table
# with NO rows (a genome with no profiled sample at this timepoint) has
# exactly the header of a table with rows, and every reader finds its columns.
# ``process_mag_files`` fills one value per name; a test keeps the two in step.
QC_SAMPLE_COLUMNS = (
    "sample_id",
    "MAG_ID",
    "file_path",
    "group",
    "subjectID",
    "replicate",
    "genome_size",
    "breadth",
    "breadth_threshold_passed",
    "breadth_fail_reason",
    "average_coverage",
    "average_coverage_including_zeros",
    "coverage_threshold_passed",
    "coverage_fail_reason",
    "median_coverage",
    "median_coverage_including_zeros",
    "coverage_std",
    "coverage_std_including_zeros",
    "length_weighted_coverage",
)


def qc_table_columns(data_type):
    """The header of a per-genome QC table for this data type.

    ``QC_SAMPLE_COLUMNS``, plus ``time`` for longitudinal data (the metadata
    of single-timepoint data carries no timepoint).

    Example: ``qc_table_columns("single")`` -> the 19 per-sample columns;
    ``qc_table_columns("longitudinal")`` -> the same 19, then ``"time"``.
    """
    columns = list(QC_SAMPLE_COLUMNS)
    if data_type == "longitudinal":
        columns.append("time")
    return columns


def init_worker(
    metadata,
    mag_sizes,
    contig_lengths,
    contig_to_mag_map,
    mag_to_contigs_map,
):
    """Initialize worker process globals shared across all sample tasks."""
    global metadata_dict
    global mag_size_dict
    global contig_length_dict
    global contig_to_mag
    global mag_to_contigs
    metadata_dict = metadata
    mag_size_dict = mag_sizes
    contig_length_dict = contig_lengths
    contig_to_mag = contig_to_mag_map
    mag_to_contigs = mag_to_contigs_map


def process_mag_files(args):
    """
    Process a MAG (Metagenome-Assembled Genome) profile file and calculate coverage statistics.

    This function reads a coverage profile for a specific MAG in a sample and determines
    whether it passes the breadth of coverage threshold and coverage depth threshold.
    The function uses global dictionaries `metadata_dict` and `mag_size_dict` to access
    sample metadata and MAG sizes.

    Parameters
    ----------
    args : tuple
        A tuple containing:
        - sample_id (str): Identifier for the sample
        - profile_fPath (str): Path to the coverage profile file
        - mag_id (str): Identifier for the MAG
        - breadth_threshold (float): Minimum required breadth of coverage (0.0-1.0)
        - coverage_threshold (float): Minimum required average coverage depth

    Returns
    -------
    dict
        A dictionary containing coverage statistics and metadata:
        - sample_id: Sample identifier
        - MAG_ID: MAG identifier
        - file_path: Path to the analyzed profile
        - group: Sample group from metadata
        - subjectID: Subject identifier from metadata
        - replicate: Replicate information from metadata
        - time: Time point (if available in metadata)
        - genome_size: Size of the MAG in base pairs
        - breadth: Calculated breadth of coverage (proportion of genome covered)
        - breadth_threshold_passed: Boolean indicating if breadth threshold was met
        - breadth_fail_reason: Explanation if threshold not met
        - average_coverage: Mean depth over COVERED positions (MIDAS2's mean genome
          coverage; the per-site depth filter's reference, computed once here)
        - average_coverage_including_zeros: Summed depth / genome size (absent
          positions as zero); the number the coverage threshold tests
        - coverage_threshold_passed: Boolean indicating if coverage threshold was met
        - coverage_fail_reason: Explanation if threshold not met
        - median_coverage: Median of per-position coverage (A+C+G+T)
        - median_coverage_including_zeros: Median coverage including zeros for absent positions
        - coverage_std: Standard deviation of per-position coverage (population, ddof=0)
        - coverage_std_including_zeros: Population std over entire genome length (absent positions as zero)
        - length_weighted_coverage: Length-weighted mean: Σ(mean_cov_contig * contig_length)/Σ(contig_length)

    Raises
    ------
    ValueError
        If more than 2 unique groups are found in the valid samples.

    Notes
    -----
    Breadth of coverage is calculated as the number of positions with at least 1x coverage
    divided by the genome size.
    """
    sample_id, profile_fPath, mag_id, breadth_threshold, coverage_threshold = args

    # Note: QC metrics (breadth, coverage) are group-independent.
    # Group-pair specific eligibility is handled downstream by eligibility_table.py.

    # Build a dictionary to store coverage stats
    # Columns are organized by category for better readability
    result = {
        # Sample metadata
        "sample_id": sample_id,
        "MAG_ID": mag_id,
        "file_path": profile_fPath,
        "group": metadata_dict[sample_id]["group"],
        "subjectID": metadata_dict[sample_id]["subjectID"],
        "replicate": metadata_dict[sample_id]["replicate"],
        "genome_size": mag_size_dict.get(mag_id, 0),
        # Breadth metrics (grouped together)
        "breadth": None,
        "breadth_threshold_passed": False,
        "breadth_fail_reason": "",
        # Coverage metrics (grouped together).  Two numbers, two jobs
        # (2026-10-15): average_coverage = covered-sites mean, the site-depth
        # filter's reference; average_coverage_including_zeros = genome-size
        # mean, the QC threshold's number.
        "average_coverage": None,
        "average_coverage_including_zeros": None,
        "coverage_threshold_passed": False,
        "coverage_fail_reason": "",
        # Additional coverage statistics
        "median_coverage": None,
        "median_coverage_including_zeros": None,
        "coverage_std": None,
        "coverage_std_including_zeros": None,
        "length_weighted_coverage": None,
    }

    # Add time if available in metadata
    if "time" in metadata_dict[sample_id]:
        result["time"] = metadata_dict[sample_id]["time"]

    # Check if MAG size is known
    mag_size = mag_size_dict.get(mag_id)
    if mag_size is None or mag_size <= 0:
        msg = f"MAG size for {mag_id} not found or invalid."
        logger.error(f"{msg} sample={sample_id}, profile_fPath={profile_fPath}")
        result["breadth_fail_reason"] = msg
        return result

    # Load and validate profile (filter to MAG's contigs).  total_coverage is
    # A+C+G+T by construction (profile_mags.build_profile_row), the same number
    # the site-depth filter later divides by, so it is read as is.
    try:
        df = load_and_validate_profile(
            profile_fPath,
            mag_id,
            sample_id,
            contig_to_mag,
            usecols=["contig", "total_coverage"],
            mag_contigs=mag_to_contigs.get(mag_id),
        )
    except ValueError as e:
        # Profile empty after filtering to MAG contigs
        result["breadth_fail_reason"] = str(e)
        return result

    # Check for zero coverage values in the profile (should not exist: the
    # writer skips positions with no A/C/G/T read)
    zero_coverage_count = (df["total_coverage"] == 0).sum()
    if zero_coverage_count > 0:
        logger.warning(
            f"Profile contains {zero_coverage_count} positions with zero coverage for MAG {mag_id}. "
            f"These should typically be absent from the profile. file={profile_fPath} sample={sample_id}"
        )

    # Use genome size as denominator for all metrics
    denom = mag_size

    # Calculate breadth metrics (genome-wide only)
    breadth_metrics = calculate_breadth_metrics(
        df, denom, mag_size, using_filter=False, positions_with_coverage_genome=None
    )
    # Only update with non-None values (filters out breadth_genome)
    result.update({k: v for k, v in breadth_metrics.items() if v is not None})

    # Calculate coverage metrics (genome-wide only)
    coverage_metrics = calculate_coverage_metrics(
        df, denom, mag_size, using_filter=False, total_coverage_sum_genome=None
    )
    # Only update with non-None values (filters out average_coverage_genome)
    result.update({k: v for k, v in coverage_metrics.items() if v is not None})

    # The including-zeros mean is the one the whole-genome std formula and the
    # QC threshold are defined on (the covered-sites mean would pass any
    # threshold for a genome with a few deep positions).
    avg_cov_including_zeros = result["average_coverage_including_zeros"]

    # Calculate median and standard deviation metrics
    median_std_metrics = calculate_median_and_std_metrics(
        df, denom, avg_cov_including_zeros
    )
    result.update(median_std_metrics)

    # Calculate length-weighted coverage
    length_weighted_cov = calculate_length_weighted_coverage(
        df, contig_length_dict, mag_id, sample_id, mag_size
    )
    result["length_weighted_coverage"] = length_weighted_cov

    # Check breadth threshold
    breadth_passed, breadth_fail_reason = check_breadth_threshold(
        result["breadth"],
        breadth_genome=None,
        threshold=breadth_threshold,
        using_filter=False,
        sample_id=sample_id,
        mag_id=mag_id,
    )
    result["breadth_threshold_passed"] = breadth_passed
    result["breadth_fail_reason"] = breadth_fail_reason

    # Check coverage threshold (cascades from breadth)
    coverage_passed, coverage_fail_reason = check_coverage_threshold(
        avg_cov_including_zeros,
        avg_cov_genome=None,
        threshold=coverage_threshold,
        breadth_passed=breadth_passed,
        using_filter=False,
        sample_id=sample_id,
        mag_id=mag_id,
    )
    result["coverage_threshold_passed"] = coverage_passed
    result["coverage_fail_reason"] = coverage_fail_reason

    return result


def _aggregate_mag(mag_id, results_list, output_dir, data_type):
    """Write one genome's QC table and build its three summaries.

    Runs in the main process after all samples for this MAG have been
    processed by the pool.  Pure pandas operations on a small DataFrame.

    The table holds the per-sample results as they are, one row per sample, in
    the fixed column order of ``qc_table_columns``.  A genome with no sample
    (``results_list`` empty: its metadata file was header-only) gets the header
    and no rows, so "no sample of this genome here" is a file on disk and a
    missing table can only mean a failed step.

    Example: two samples of ``MAG_A`` -> ``MAG_A_QC.tsv`` with two rows and the
    summaries of the samples that passed; no sample of ``MAG_B`` ->
    ``MAG_B_QC.tsv`` with the header only and ``(None, None, None)``.

    Returns ``(group_summary, group_time_summary, overall_summary)``; all None
    when the genome has no sample or none passed QC.
    """
    columns = qc_table_columns(data_type)
    out_file = Path(output_dir) / f"{mag_id}_QC.tsv"

    if not results_list:
        pd.DataFrame(columns=columns).to_csv(out_file, sep="\t", index=False)
        logger.info(
            f"MAG {mag_id}: no sample in its metadata file; "
            f"header-only QC table written: {out_file}"
        )
        return None, None, None

    df_results = pd.DataFrame(results_list)
    # The sample workers must return exactly the table's columns: a silent
    # extra or missing one would give tables with rows a different header from
    # the header-only ones.
    if set(df_results.columns) != set(columns):
        raise ValueError(
            f"MAG {mag_id}: QC results have columns {sorted(df_results.columns)}, "
            f"expected {sorted(columns)}"
        )
    df_results = df_results[columns]

    df_results.to_csv(out_file, sep="\t", index=False)
    logger.debug(f"Saved QC report for {mag_id} to {out_file}")

    if not df_results["coverage_threshold_passed"].any():
        logger.warning(f"MAG {mag_id}: no samples passed QC; skipping summaries.")
        return None, None, None

    # Prepare summaries
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
    setup_logging()

    parser = argparse.ArgumentParser(
        description=("Script to perform genome-wide quality control on MAGs."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rootDir",
        required=True,
        type=Path,
        help="Path to the root directory containing metadata files.",
    )
    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="FASTA file with contigs (for MAG size)",
    )
    parser.add_argument(
        "--breadth_threshold",
        type=float,
        default=0.1,
        help="Minimum breadth coverage required to pass.",
    )
    parser.add_argument(
        "--coverage_threshold",
        type=float,
        default=1,
        help="Minimum average coverage depth required to pass.",
    )
    parser.add_argument(
        "--cpus", type=int, default=cpu_count(), help="Number of processors to use."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=Path,
        help="Directory to write output table.",
    )

    parser.add_argument(
        "--data_type",
        help="Is the data from a single timepoint or from a time series (longitudinal)",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )
    parser.add_argument(
        "--mag_mapping_file",
        help="Path to tab-separated file mapping contig names to MAG IDs. "
        "Must have columns 'contig_id' and 'mag_id'.",
        type=Path,
        required=True,
        metavar="filepath",
    )

    args = parser.parse_args()

    # Note: For position-based QC analysis, use the separate positions_qc.py script
    # located in alleleflux/scripts/accessory/positions_qc.py

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Using breadth threshold: {args.breadth_threshold}")
    logger.info(f"Using coverage threshold: {args.coverage_threshold}")

    mag_size_dict = calculate_mag_sizes(args.fasta, args.mag_mapping_file)
    # Build contig length index once (re-used by workers)
    contig_length_dict = build_contig_length_index(args.fasta, args.mag_mapping_file)
    # Load mapping for MAG membership validation/subsetting
    contig_to_mag_global = load_mag_mapping(args.mag_mapping_file)

    metadata_files = list(args.rootDir.glob("*_metadata.tsv"))
    if not metadata_files:
        raise ValueError("No *_metadata.tsv files found in input directory")

    # ── Build a flat list of (sample) tasks and a per-MAG metadata index ──
    # load_mag_metadata_file returns:
    #   metadata_dict: {sample_id -> {group, subjectID, replicate, [time]}}
    #   sample_tuples: [(sample_id, profile_path, mag_id, breadth, cov), ...]
    #
    # We merge ALL per-MAG metadata_dicts into one combined dict (sample_ids
    # are globally unique) and build a flat task list so the Pool dispatches
    # at the sample level, giving perfect work-stealing across MAGs of
    # different sizes.
    combined_metadata = {}  # sample_id -> metadata row
    mag_sample_ids = {}  # mag_id -> [sample_id, ...] (for post-pool grouping)
    all_sample_tasks = []  # flat list fed to pool.imap_unordered

    for metadata_file in metadata_files:
        mag_id = metadata_file.stem.split("_metadata")[0]
        metadata_dict, sample_tuples = load_mag_metadata_file(
            metadata_file,
            mag_id,
            args.breadth_threshold,
            args.coverage_threshold,
            args.data_type,
        )
        # A header-only metadata file (the genome is profiled elsewhere in the
        # run but has no sample at this timepoint) gives no task; the genome
        # keeps its entry so that it still gets a header-only QC table below.
        combined_metadata.update(metadata_dict)
        mag_sample_ids[mag_id] = [t[0] for t in sample_tuples]
        all_sample_tasks.extend(sample_tuples)

    # One genome without a sample is normal; NO sample in the whole folder is
    # not: there is nothing to read, and a folder of header-only QC tables
    # would let the run carry on with no data.
    if not all_sample_tasks:
        raise ValueError(
            f"no sample in any of the {len(metadata_files)} metadata file(s) of "
            f"{args.rootDir}: every file is header-only, so there is nothing to QC"
        )

    n_workers = min(args.cpus, len(all_sample_tasks))
    n_mags = len(mag_sample_ids)
    logger.info(
        f"Processing {len(all_sample_tasks)} samples across {n_mags} MAGs "
        f"with {n_workers} parallel workers."
    )

    # Invert contig_to_mag once so load_and_validate_profile can do an O(1)
    # membership test per unique contig instead of an O(all_contigs) map per row.
    from collections import defaultdict

    mag_to_contigs_global = defaultdict(set)
    for contig, mid in contig_to_mag_global.items():
        mag_to_contigs_global[mid].add(contig)

    # Silence verbose utility loggers before forking workers.
    logging.getLogger("alleleflux.scripts.utilities").setLevel(logging.WARNING)

    # ── Run sample-level pool ──
    # The large dicts (combined_metadata, mag_size_dict, contig_length_dict,
    # contig_to_mag_global, mag_to_contigs_global) are sent ONCE per worker via
    # the initializer, not once per task.
    with Pool(
        processes=n_workers,
        initializer=init_worker,
        initargs=(
            combined_metadata,
            mag_size_dict,
            contig_length_dict,
            contig_to_mag_global,
            mag_to_contigs_global,
        ),
    ) as pool:
        with logging_redirect_tqdm():
            sample_results = list(
                tqdm(
                    pool.imap_unordered(
                        process_mag_files, all_sample_tasks, chunksize=4
                    ),
                    total=len(all_sample_tasks),
                    desc="Processing samples",
                    unit="sample",
                )
            )

    # ── Group results by MAG and run per-MAG aggregation ──
    # This is pure pandas on small DataFrames — fast in the main process.
    results_by_mag = {mag_id: [] for mag_id in mag_sample_ids}
    for result in sample_results:
        results_by_mag[result["MAG_ID"]].append(result)

    all_group = []
    all_group_time = []
    all_overall = []

    for mag_id in tqdm(mag_sample_ids, desc="Aggregating MAGs", unit="MAG"):
        group_df, group_time_df, overall_df = _aggregate_mag(
            mag_id,
            results_by_mag[mag_id],
            args.output_dir,
            args.data_type,
        )
        if group_df is not None:
            all_group.append(group_df)
        if group_time_df is not None:
            all_group_time.append(group_time_df)
        if overall_df is not None:
            all_overall.append(overall_df)

    # Write combined summaries
    summaries = [
        (all_group, "ALL_MAGs_QC_group_summary.tsv"),
        (all_group_time, "ALL_MAGs_QC_group_time_summary.tsv"),
        (all_overall, "ALL_MAGs_QC_overall_summary.tsv"),
    ]
    for dfs, filename in summaries:
        if dfs:
            pd.concat(dfs, ignore_index=True).to_csv(
                args.output_dir / filename, sep="\t", index=False
            )

    logger.info("QC analysis completed.")


if __name__ == "__main__":
    main()
