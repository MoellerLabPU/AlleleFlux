#!/usr/bin/env python3

import argparse
import logging
import sys
from multiprocessing import Pool, cpu_count
from pathlib import Path

import numpy as np
import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.utilities import (
    build_contig_length_index,
    calculate_mag_sizes,
    load_mag_mapping,
    load_mag_metadata_file,
)

logger = logging.getLogger(__name__)


def load_positions_file(positions_file_path):
    """
    Load positions filter file and return MAG -> set of (contig, position) tuples.

    Parameters:
        positions_file_path (Path or str): Path to TSV file with columns: MAG, contig, position

    Returns:
        dict: Mapping of MAG ID -> set of (contig, position) tuples

    Raises:
        ValueError: If duplicate positions are found for a MAG
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


def init_worker(
    metadata,
    mag_sizes,
    contig_lengths,
    contig_to_mag_map,
    positions_filter_map=None,
    positions_denominator="positions",
):
    """
    Initialize worker process with metadata and MAG sizes.

    This function sets up global dictionaries for metadata and MAG sizes
    that can be accessed by worker processes.

    Parameters:
        metadata (dict): A dictionary containing metadata information.
        mag_sizes (dict): A dictionary containing sizes of MAGs (Metagenome-Assembled Genomes).
        positions_filter_map (dict): Optional mapping MAG -> set of (contig, position) tuples.
        positions_denominator (str): Either 'positions' or 'genome' - determines denominator for metrics.

    Returns:
        None
    """
    global metadata_dict
    global mag_size_dict
    metadata_dict = metadata
    mag_size_dict = mag_sizes
    # Store contig lengths for weighted coverage
    global contig_length_dict
    contig_length_dict = contig_lengths
    # Mapping: contig_id -> mag_id (for validation/subsetting)
    global contig_to_mag
    contig_to_mag = contig_to_mag_map
    # Optional: mapping MAG -> DataFrame of allowed (contig, position)
    global positions_filter
    positions_filter = positions_filter_map
    # Denominator mode: 'positions' or 'genome'
    global positions_den_mode
    positions_den_mode = positions_denominator


def process_mag_files(args):
    """
    Process a MAG (Metagenome-Assembled Genome) profile file and calculate coverage statistics.

    This function reads a coverage profile for a specific MAG in a sample and determines
    whether it passes the breadth of coverage threshold. The function uses global dictionaries
    `metadata_dict` and `mag_size_dict` to access sample metadata and MAG sizes.

    Parameters
    ----------
    args : tuple
        A tuple containing:
        - sample_id (str): Identifier for the sample
        - profile_fPath (str): Path to the coverage profile file
        - mag_id (str): Identifier for the MAG
        - breadth_threshold (float): Minimum required breadth of coverage (0.0-1.0)

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
        - average_coverage: Average genome coverage depth (sum of all total_coverage values divided by genome size)
        - median_coverage: Median of per-position total_coverage
        - length_weighted_coverage: Length-weighted mean of per-contig mean coverage: Σ(mean_cov_contig * contig_length)/Σ(contig_length)
        - coverage_std: Standard deviation of per-position total_coverage (population, ddof=0)
        - coverage_std_including_zeros: Population std over entire genome length (absent positions treated as zero)
        - breadth_threshold_passed: Boolean indicating if breadth threshold was met
        - breadth_fail_reason: Explanation if threshold not met

    Raises
    ------
    ValueError
        If more than 2 unique groups are found in the valid samples.

    Notes
    -----
    Breadth of coverage is calculated as the number of positions with at least
    1x coverage divided by the total genome size.
    """
    sample_id, profile_fPath, mag_id, breadth_threshold = args

    # Check if there are more than 2 unique groups in the metadata
    unique_groups = {meta["group"] for meta in metadata_dict.values()}
    if len(unique_groups) > 2:
        msg = f"More than 2 unique groups found: {unique_groups}"
        raise ValueError(msg)

    # Check if we'll be using a positions filter (determines which fields are added)
    using_position_filter = (
        positions_filter is not None
        and positions_filter.get(mag_id) is not None
        and len(positions_filter.get(mag_id, set())) > 0
    )

    # Build a dictionary to store coverage stats
    result = {
        "sample_id": sample_id,
        "MAG_ID": mag_id,
        "file_path": profile_fPath,
        "group": metadata_dict[sample_id]["group"],
        "subjectID": metadata_dict[sample_id]["subjectID"],
        "replicate": metadata_dict[sample_id]["replicate"],
        "genome_size": mag_size_dict.get(mag_id, 0),
        "breadth": None,
        "average_coverage": None,
        "median_coverage": None,
        "median_coverage_including_zeros": None,
        "coverage_std": None,
        "coverage_std_including_zeros": None,
    }

    # Add breadth_genome field when dual breadth mode is active
    # (positions file provided AND positions_denominator="positions")
    if using_position_filter and positions_den_mode == "positions":
        result["breadth_genome"] = None

    # Only include length_weighted_coverage when positions filter is NOT active
    # (this metric is not meaningful for sparse position sets)
    if not using_position_filter:
        result["length_weighted_coverage"] = None

    if using_position_filter:
        logger.debug(
            f"Positions filter detected for MAG {mag_id} with {len(positions_filter.get(mag_id, set()))} positions"
        )

    # Add breadth threshold fields in two scenarios:
    # 1. When positions filter is NOT active (standard case)
    # 2. When positions filter is active AND positions_denominator="positions" (dual breadth mode)
    if (
        positions_filter is None
        or not using_position_filter
        or (using_position_filter and positions_den_mode == "positions")
    ):
        result["breadth_threshold_passed"] = False
        result["breadth_fail_reason"] = ""
    # Add time if available in metadata
    if "time" in metadata_dict[sample_id]:
        result["time"] = metadata_dict[sample_id]["time"]

    # Check if MAG size is known
    mag_size = mag_size_dict.get(mag_id)
    if mag_size is None or mag_size <= 0:
        msg = f"MAG size for {mag_id} not found or invalid."
        logger.error(f"{msg} sample={sample_id}, profile_fPath={profile_fPath}")
        if positions_filter is None or not using_position_filter:
            result["breadth_fail_reason"] = msg
        return result

    df = pd.read_csv(profile_fPath, sep="\t", dtype={"gene_id": str})

    # Ensure the profile only includes contigs for this MAG; log if rows are dropped
    if "contig" in df.columns:
        before_rows = len(df)
        # Map contigs to mag; unknown contigs map to None and will be excluded
        mag_lookup = df["contig"].map(contig_to_mag)
        mask = mag_lookup == mag_id
        dropped = int(before_rows - mask.sum())
        if dropped > 0:
            logger.warning(
                f"Profile contains {dropped} rows not belonging to MAG {mag_id}; filtering them out. "
                f"file={profile_fPath} sample={sample_id}"
            )
        df = df.loc[mask].copy()
        if df.empty:
            msg = (
                f"No rows remain after filtering profile to MAG {mag_id}. "
                f"file={profile_fPath}"
            )
            logger.error(msg)
            if positions_filter is None or not using_position_filter:
                result["breadth_fail_reason"] = msg
            return result

    # In dual breadth mode, capture genome-wide positions with coverage BEFORE filtering
    positions_with_coverage_genome = None
    if using_position_filter and positions_den_mode == "positions":
        positions_with_coverage_genome = df[df["total_coverage"] >= 1].shape[0]

    # Optionally filter to positions-of-interest if provided for this MAG
    positions_universe_n = None
    if "contig" in df.columns and "position" in df.columns and using_position_filter:
        wanted = positions_filter.get(mag_id)
        if wanted is not None and len(wanted) > 0:
            positions_universe_n = len(wanted)
            # logger.info(
            #     f"Positions filter active for {mag_id} (mode={positions_den_mode}, N={positions_universe_n})"
            # )
            # Filter df to only rows present in positions set (fast set membership)
            df["contig"] = df["contig"].astype(str)
            df["position"] = pd.to_numeric(df["position"]).astype("Int64")
            # before_rows_pos = len(df)
            # Create boolean mask using set membership (much faster than merge)
            mask = df.apply(
                lambda row: (row["contig"], int(row["position"])) in wanted, axis=1
            )
            df = df.loc[mask].copy()
            # dropped_pos = before_rows_pos - len(df)
            # if dropped_pos > 0:
            #     logger.info(
            #         f"Filtered out {dropped_pos} positions not in provided list for MAG {mag_id} (sample {sample_id})."
            #     )
            if df.empty:
                msg = f"No matching (contig, position) rows for MAG {mag_id} in sample {sample_id} after applying positions filter."
                logger.info(msg)
                result["positions_considered"] = positions_universe_n or 0
                # Note: breadth_threshold_passed and breadth_fail_reason fields don't exist when using filter
                return result
    elif using_position_filter:
        # Position filter was intended but profile doesn't have required columns
        logger.warning(
            f"Positions filter specified for {mag_id} but profile file lacks 'contig' and/or 'position' columns. "
            f"Available columns: {df.columns.tolist()}. Sample: {sample_id}"
        )

    # Check for zero coverage values in the profile (should not exist)
    zero_coverage_count = (df["total_coverage"] == 0).sum()
    if zero_coverage_count > 0:
        logger.warning(
            f"Profile contains {zero_coverage_count} positions with zero coverage for MAG {mag_id}. "
            f"These should typically be absent from the profile. file={profile_fPath} sample={sample_id}"
        )

    # Compute metrics denominator: genome_size (default) or positions universe if filtered
    if using_position_filter and positions_universe_n:
        denom = positions_universe_n if positions_den_mode == "positions" else mag_size
        # Only include positions_considered when a positions file is provided
        result["positions_considered"] = positions_universe_n
    else:
        denom = mag_size

    # Compute breadth coverage
    positions_with_coverage = df[df["total_coverage"] >= 1].shape[0]
    breadth = positions_with_coverage / denom if denom > 0 else np.nan
    result["breadth"] = breadth

    # Compute breadth_genome when dual breadth mode is active
    # (positions file provided AND positions_denominator="positions")
    if using_position_filter and positions_den_mode == "positions":
        breadth_genome = (
            positions_with_coverage_genome / mag_size if mag_size > 0 else np.nan
        )
        result["breadth_genome"] = breadth_genome

    # Average coverage depth = sum of coverage across all rows divided by denominator
    total_coverage_sum = df["total_coverage"].sum()
    avg_cov = total_coverage_sum / denom if denom > 0 else np.nan
    result["average_coverage"] = avg_cov

    # Median coverage (per-position) among observed covered bases only (exclude zeros if any slipped in)
    pos_cov = df.loc[df["total_coverage"] > 0, "total_coverage"]
    result["median_coverage"] = float(pos_cov.median()) if not pos_cov.empty else np.nan

    # Median coverage including zeros for absent positions across entire genome
    observed_positions = len(df)
    absent_positions = int(denom - observed_positions) if denom is not None else 0
    if absent_positions >= 0:
        # Create array with observed coverage values plus zeros for absent positions
        all_coverage = np.concatenate(
            [df["total_coverage"].values, np.zeros(absent_positions)]
        )
        result["median_coverage_including_zeros"] = float(np.median(all_coverage))
    else:
        # This shouldn't happen, but handle gracefully
        logger.warning(
            f"Profile has more positions ({observed_positions}) than denominator ({denom}) for {mag_id} sample {sample_id}"
        )
        result["median_coverage_including_zeros"] = result["median_coverage"]

    # Standard deviation across covered positions only (population std, ddof=0)
    result["coverage_std"] = float(pos_cov.std(ddof=0)) if not pos_cov.empty else np.nan

    # Population std including zeros across entire genome length (absent positions treated as 0)
    if denom > 0 and np.isfinite(avg_cov):
        ex2 = float((df["total_coverage"] ** 2).sum()) / denom
        var = max(ex2 - (avg_cov * avg_cov), 0.0)
        result["coverage_std_including_zeros"] = float(np.sqrt(var))
    else:
        result["coverage_std_including_zeros"] = np.nan

    # --- Length-Weighted Coverage Metrics (contig-level) ---
    # Only compute when positions filter is NOT active (not meaningful for sparse positions)
    if not using_position_filter:
        # Expect a column 'contig' with per-position coverage; aggregate to per-contig mean first.
        # Sum coverage per contig; later divide by contig length to get mean over all bases
        contig_means = (
            df.groupby("contig", as_index=False)["total_coverage"]
            .sum()
            .rename(columns={"total_coverage": "coverage_sum"})
        )
        if not contig_means.empty:
            # Map lengths (no need to declare global again; read-only lookup)
            contig_means["length"] = contig_means["contig"].map(contig_length_dict)
            before = contig_means.shape[0]
            contig_means = contig_means[
                pd.to_numeric(contig_means["length"], errors="coerce") > 0
            ]
            dropped = before - contig_means.shape[0]
            if dropped > 0:
                logger.warning(
                    f"Dropped {dropped} contigs without valid lengths for {mag_id} (sample {sample_id})."
                )
            if not contig_means.empty:
                # Compute per-contig mean as sum(coverage)/contig_length to include implicit zeros
                contig_means["mean_coverage"] = (
                    contig_means["coverage_sum"] / contig_means["length"]
                )
                total_len = contig_means["length"].sum()
                if total_len > 0:
                    weighted_sum = (
                        contig_means["mean_coverage"] * contig_means["length"]
                    ).sum()
                    result["length_weighted_coverage"] = weighted_sum / total_len

    # Perform breadth threshold check in two scenarios:
    # 1. When positions filter is NOT active (standard case) - check against breadth
    # 2. When positions filter is active AND positions_denominator="positions" (dual breadth mode) - check against breadth_genome
    if positions_filter is None or not using_position_filter:
        # Standard case: check breadth (genome-based or full genome without filter)
        if breadth < breadth_threshold:
            msg = f"Breadth {breadth:.2%} < threshold {breadth_threshold:.2%}"
            logger.info(f"{msg} - sample={sample_id}, mag={mag_id}")
            result["breadth_fail_reason"] = msg
        else:
            result["breadth_threshold_passed"] = True
    elif using_position_filter and positions_den_mode == "positions":
        # Dual breadth mode: check breadth_genome against threshold
        breadth_genome_val = result.get("breadth_genome", np.nan)
        if breadth_genome_val < breadth_threshold:
            msg = f"Breadth (genome) {breadth_genome_val:.2%} < threshold {breadth_threshold:.2%}"
            logger.info(f"{msg} - sample={sample_id}, mag={mag_id}")
            result["breadth_fail_reason"] = msg
        else:
            result["breadth_threshold_passed"] = True

    return result


def check_timepoints(df, data_type):
    """
    Check if MAGs (Metagenome-Assembled Genomes) are present in two timepoints for each subject.
    This function verifies that subjects have data from exactly two timepoints for longitudinal analysis.
    For single datatype, it skips the timepoint checks and sets two_timepoints_passed to match breadth_threshold_passed.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing sample data with 'subjectID', 'time', and 'breadth_threshold_passed' columns
    data_type : str
        Type of data analysis, either "single" or assumed to be longitudinal

    Returns
    -------
    pandas.DataFrame
        The input dataframe with an additional 'two_timepoints_passed' column indicating samples
        that come from subjects with exactly 2 timepoints

    Raises
    ------
    ValueError
        If more than 2 unique timepoints are found in the passed samples

    Notes
    -----
    - For longitudinal analysis, the function ensures there are exactly 2 unique timepoints globally
    - A sample is marked as passing if its subject has data at both timepoints and the sample has passed
      the breadth threshold
    - Logs information about subjects with valid timepoints
    """
    # For single datatype, skip timepoint checks
    if data_type == "single":
        # Set two_timepoints_passed to match breadth_threshold_passed for single data type
        df["two_timepoints_passed"] = df["breadth_threshold_passed"]
        return df

    # Initialize columns
    df["two_timepoints_passed"] = False

    # Get samples that passed coverage threshold
    passed_samples = df[df["breadth_threshold_passed"]]

    # Check that there are exactly 2 unique timepoints globally.
    unique_times = passed_samples["time"].unique()
    if len(unique_times) > 2:
        raise ValueError(f"More than 2 unique timepoints found: {unique_times}")

    if passed_samples.empty:
        logger.warning("No samples passed coverage threshold")
        return df

    # Find subjects with MAG in exactly 2 timepoints
    valid_subjects = (
        passed_samples.groupby("subjectID")["time"]
        .nunique()
        .where(lambda x: x == 2)
        .dropna()
        .index
    )

    # Mark passed samples from valid subjects
    valid_mask = passed_samples["subjectID"].isin(valid_subjects)
    df.loc[passed_samples[valid_mask].index, "two_timepoints_passed"] = True

    # Add debug info
    if not valid_subjects.empty:
        logger.info(f"Subjects with 2 valid timepoints: {len(valid_subjects)}")
    else:
        logger.warning("No subjects have valid MAG in 2 timepoints")

    return df


def add_subject_count_per_group(df):
    """
    Add columns for subject count and replicate count per group to the DataFrame.

    This function calculates the number of unique subjects and replicates for each group
    in the dataset, but only considers samples that have passed timepoint validation
    (where 'two_timepoints_passed' is True). It then adds these counts as new columns
    to the original DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least the following columns:
        - 'two_timepoints_passed': Boolean indicating if the sample has passed timepoint validation
        - 'group': Group identifier
        - 'subjectID': Subject identifier
        - 'replicate': Replicate identifier

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with two new columns added:
        - 'subjects_per_group': Number of unique subjects in the group (NaN for invalid samples)
        - 'replicates_per_group': Number of unique replicates in the group (NaN for invalid samples)

    Notes
    -----
    If no valid samples are found, both new columns will be set to NaN for all rows.
    The function logs a summary of subject and replicate counts per group for valid samples.
    """
    # Get samples that passed timepoint validation
    valid_subjects = df[df["two_timepoints_passed"]]

    if valid_subjects.empty:
        logger.warning("No valid subjects found")
        # Initialize the columns with NaN values
        df["subjects_per_group"] = np.nan
        df["replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "subjects_per_group"] = 0
        df.loc[df["two_timepoints_passed"], "replicates_per_group"] = 0
        return df

    # Count unique subjects per group in valid samples
    group_counts = (
        valid_subjects.groupby("group")["subjectID"]
        .nunique()  # Count unique subjects per group
        .reset_index(name="subjects_per_group")
    )

    # Count replicates per group in valid samples
    replicate_counts = (
        valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="replicates_per_group")
    )

    # Merge both counts into one counts DataFrame keyed by group.
    counts = group_counts.merge(replicate_counts, on="group", how="outer")

    # Merge the counts back into the original DataFrame based on group
    df = df.merge(counts, on="group", how="left")

    # For rows where two_timepoints_passed is False, set the counts to NaN.
    df.loc[
        ~df["two_timepoints_passed"], ["subjects_per_group", "replicates_per_group"]
    ] = np.nan

    # Log summary
    logger.info("Subject and replicate counts per group (valid samples):")
    logger.info(f"\n{counts}")
    return df


def count_paired_replicates(df):
    """
    Count the number of paired replicates per group in the dataset.

    This function identifies valid subjects (those that passed both timepoints) and counts
    replicates that appear in both groups. It ensures exactly two groups exist in the valid
    samples dataset.

    A "paired replicate" means the same replicate has samples in both experimental groups.
    This helps identify balanced experimental designs where the same experimental units
    are measured across different conditions.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing the following columns:
        - 'two_timepoints_passed': Boolean indicating if the sample has two valid timepoints
        - 'group': The experimental group identifier
        - 'replicate': The replicate identifier for the sample

    Returns
    -------
    pandas.DataFrame
        The original DataFrame with an additional column:
        - 'paired_replicates_per_group': Number of replicates per group that have samples in both groups.
          NaN for samples that failed timepoint check or whose replicate isn't paired.
          Zero (0) for samples that passed timepoint check but have no paired replicates.


    Notes
    -----
    - If there are no valid subjects at all, the function returns the DataFrame with
      'paired_replicates_per_group' column set to NaN for all rows.
    - If there are valid subjects but no paired replicates, the function sets
      'paired_replicates_per_group' to 0 for valid subjects (those that passed timepoint check).
    - For samples that failed the timepoint check, the column is always set to NaN.
    """
    valid_subjects = df[df["two_timepoints_passed"]]

    if valid_subjects.empty:
        logger.warning("No valid subjects found")
        # Initialize the column with NaN values
        df["paired_replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "paired_replicates_per_group"] = 0
        return df

    rep_group_count = (
        valid_subjects.groupby("replicate")["group"]
        .nunique()
        .reset_index(name="group_count")
    )
    paired_reps = rep_group_count[rep_group_count["group_count"] == 2]["replicate"]

    # Check if there are any paired replicates
    if paired_reps.empty:
        logger.warning("No paired replicates found across groups")
        # Initialize the column with NaN values
        df["paired_replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "paired_replicates_per_group"] = 0
        return df

    paired_valid_subjects = valid_subjects[
        valid_subjects["replicate"].isin(paired_reps)
    ]

    # Count the number of unique paired replicates per group
    paired_counts = (
        paired_valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="paired_replicates_per_group")
    )

    # Merge the paired replicate counts back into the original DataFrame based on 'group'
    df = df.merge(paired_counts, on="group", how="left")

    # Set counts to NaN where:
    # - Samples failed timepoint check OR
    # - Samples passed timepoint check but replicate isn't paired
    valid_mask = df["two_timepoints_passed"]
    paired_mask = df["replicate"].isin(paired_reps)
    df["paired_replicates_per_group"] = np.where(
        valid_mask & paired_mask, df["paired_replicates_per_group"], np.nan
    )

    # Log summary
    logger.info("Replicates paired per group")
    logger.info(f"\n{paired_counts}")
    return df


def _aggregate_mag_stats(
    df_results: pd.DataFrame,
    mag_id: str,
    group_cols: list,
    overall: bool = False,
):
    """
    Generic aggregator for MAG QC metrics.

    Parameters
    ----------
    df_results : pd.DataFrame
        Sample-level QC results for a single MAG.
    mag_id : str
        MAG identifier.
    group_cols : list
        Column names to group by. Empty list means collapse across all samples.
    output_filename : str
        Name of TSV to write.
    output_dir : str
        Destination directory.
    overall : bool
        If True, produce single-row overall summary with *_mean suffix; else per-group means.
    """
    if df_results.empty:
        logger.warning(f"No results to summarize for {mag_id}")
        return None

    # If breadth_threshold_passed is absent (positions mode), consider all rows as passing
    if "breadth_threshold_passed" in df_results.columns:
        passing = df_results[df_results["breadth_threshold_passed"]]
        # logger.info("breadth_threshold_passed column found.")
        # logger.info(f"Number of passing samples for {mag_id}: {passing.shape[0]}")
    else:
        # logger.info("breadth_threshold_passed column not found; using all samples.")
        passing = df_results.copy()
    if passing.empty:
        logger.warning(
            f"No passing samples for {mag_id}; summarizing all samples instead."
        )
        passing = df_results.copy()

    metric_cols = [
        c
        for c in [
            "subjects_per_group",
            "replicates_per_group",
            "paired_replicates_per_group",
            "breadth",
            "average_coverage",
            "median_coverage",
            "median_coverage_including_zeros",
            "length_weighted_coverage",
            "coverage_std",
            "coverage_std_including_zeros",
        ]
        if c in passing.columns
    ]

    if overall:
        row = {"MAG_ID": mag_id, "num_samples": passing.shape[0]}
        for col in metric_cols:
            row[col + "_mean"] = passing[col].mean()
        out_df = pd.DataFrame([row])
    else:
        if not group_cols:
            logger.error(
                "group_cols empty but overall flag False; nothing to aggregate."
            )
            return None
        agg_dict = {c: "mean" for c in metric_cols}
        agg_dict["sample_id"] = "count"
        grouped = passing.groupby(group_cols, dropna=False).agg(agg_dict)
        out_df = grouped.rename(columns={"sample_id": "num_samples"}).reset_index()
        out_df.insert(0, "MAG_ID", mag_id)

        # Add _mean suffix to all metric columns for consistency
        rename_dict = {
            col: col + "_mean" for col in metric_cols if col in out_df.columns
        }
        out_df = out_df.rename(columns=rename_dict)

    return out_df


def process_mag(args):
    (
        mag_id,
        metadata_file,
        mag_size_dict,
        contig_length_dict_global,
        contig_to_mag_map,
        positions_filter_map,
        output_dir,
        breadth_threshold,
        cpus,
        data_type,
        positions_denominator,
    ) = args
    logger.info(f"Processing MAG: {mag_id}")

    metadata_dict, sample_files_with_mag_id = load_mag_metadata_file(
        metadata_file, mag_id, breadth_threshold, data_type
    )

    if not sample_files_with_mag_id:
        raise ValueError(
            f"No samples found in the metadata file {metadata_file}. For MAG: {mag_id}. Exiting."
        )

    num_procs = min(cpus, len(sample_files_with_mag_id))
    logger.info(
        f"Processing {len(sample_files_with_mag_id)} samples for {mag_id} with {num_procs} processes."
    )

    # Load sample data in parallel
    with Pool(
        processes=num_procs,
        initializer=init_worker,
        initargs=(
            metadata_dict,
            mag_size_dict,
            contig_length_dict_global,
            contig_to_mag_map,
            positions_filter_map,
            positions_denominator,
        ),
    ) as pool:
        results_list = list(
            pool.imap_unordered(process_mag_files, sample_files_with_mag_id),
        )
    # Build results DataFrame
    df_results = pd.DataFrame(results_list)

    # When positions filter is active, skip timepoint validation and subject/replicate counting
    # (these metrics are not meaningful for subset position analysis)
    if positions_filter_map is None:
        # Check that each subject has exactly 2 timepoints (for longitudinal data only)
        df_results = check_timepoints(df_results, data_type)

        # Count the number of unique subjects and replicates per group
        df_results = add_subject_count_per_group(df_results)

        # Count the number of paired replicates per group
        df_results = count_paired_replicates(df_results)

    out_file = Path(output_dir) / f"{mag_id}_QC.tsv"
    df_results.to_csv(out_file, sep="\t", index=False)
    logger.info(f"Saved QC report for {mag_id} to {out_file}")
    # Prepare summaries (do not write here; collected in main)
    group_summary = _aggregate_mag_stats(
        df_results, mag_id, group_cols=["group"], overall=False
    )
    group_time_summary = None
    if "time" in df_results.columns:
        group_time_summary = _aggregate_mag_stats(
            df_results, mag_id, group_cols=["group", "time"], overall=False
        )
    overall_summary = _aggregate_mag_stats(
        df_results, mag_id, group_cols=[], overall=True
    )
    return group_summary, group_time_summary, overall_summary


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Script to detect which samples for a MAG pass the breadth threshold and output a pass/fail table.",
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
    parser.add_argument(
        "--positions_file",
        help=(
            "Optional TSV file with columns 'MAG', 'contig', 'position'. If provided, "
            "QC metrics are computed only for the specified positions for each MAG, and breadth/averages "
            "are calculated relative to the number of specified positions."
        ),
        type=Path,
        required=False,
        metavar="filepath",
    )
    parser.add_argument(
        "--positions_denominator",
        choices=["positions", "genome"],
        default="positions",
        help=(
            "When --positions_file is provided: use 'positions' to divide by the "
            "number of specified sites per MAG, or 'genome' to divide by genome_size."
        ),
    )

    args = parser.parse_args()

    # Validate: positions_denominator requires positions_file
    if "--positions_denominator" in sys.argv and not args.positions_file:
        raise ValueError(
            "--positions_denominator can only be used when --positions_file is provided."
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    mag_size_dict = calculate_mag_sizes(args.fasta, args.mag_mapping_file)
    # Build contig length index once (re-used by workers)
    contig_length_dict = build_contig_length_index(args.fasta, args.mag_mapping_file)
    # Load mapping for MAG membership validation/subsetting
    contig_to_mag_global = load_mag_mapping(args.mag_mapping_file)

    # Optional positions filter
    positions_by_mag = None
    if args.positions_file:
        positions_by_mag = load_positions_file(args.positions_file)
        logger.info(
            f"Loaded positions filter for {len(positions_by_mag)} MAG(s) from {args.positions_file}"
        )

    metadata_files = list(args.rootDir.glob("*_metadata.tsv"))
    if not metadata_files:
        raise ValueError("No *_metadata.tsv files found in input directory")

    # If positions file provided, only process MAGs that are in the positions file
    if positions_by_mag is not None:
        mags_in_positions = set(positions_by_mag.keys())
        metadata_files = [
            mf
            for mf in metadata_files
            if mf.stem.split("_metadata")[0] in mags_in_positions
        ]
        logger.info(
            f"Filtered to {len(metadata_files)} MAG(s) that have positions in the filter file"
        )
        if not metadata_files:
            raise ValueError(
                "No MAGs from metadata directory found in positions file. "
                "Check that MAG IDs in positions file match metadata filenames."
            )

    tasks = [
        (
            metadata_file.stem.split("_metadata")[0],  # mag_id from filename
            metadata_file,
            mag_size_dict,
            contig_length_dict,
            # pass mapping for validation
            contig_to_mag_global,
            positions_by_mag,
            args.output_dir,
            args.breadth_threshold,
            args.cpus,
            args.data_type,
            args.positions_denominator,
        )
        for metadata_file in metadata_files
    ]

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
