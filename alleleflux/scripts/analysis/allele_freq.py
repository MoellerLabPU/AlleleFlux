#!/usr/bin/env python3
"""Allele-frequency analyzer — per-combination diff/aggregate step.  [Stage 2 of 2]

Role in the two-stage pipeline
-------------------------------
This script is **Stage 2**.  It runs ONCE per (MAG, timepoint_combination,
group_combination).  Unlike Stage 1 (``allele_freq_cache.py``), which reads
raw profile TSV files, this script reads the Parquet cache files that Stage 1
already wrote and computes the combination-specific diff outputs.

    alleleflux-allele-freq
        --magID       MRGM_1079
        --cache_files allele_freq_cache/MRGM_1079_5mo_allele_frequency.parquet \\
                      allele_freq_cache/MRGM_1079_10mo_allele_frequency.parquet
        --data_type   longitudinal
        --groups      1D AL          # optional: filter cache to these groups only
        --output_dir  allele_analysis/allele_analysis_5mo_10mo-1D_AL

What changed from the original ``allele_freq.py``
--------------------------------------------------
- **Inputs**: ``--qc_file`` (single QC TSV) replaced by ``--cache_files``
  (one or two Parquet files written by ``allele_freq_cache.py``).
- **No profile reads**: all I/O is now Parquet reads — the expensive
  ``pd.read_csv`` on gzipped profile TSVs now lives entirely in Stage 1.
- **No longitudinal file written**: the old
  ``{mag}_allele_frequency_longitudinal.tsv.gz`` is no longer produced.
  Rules that previously consumed it (``lmm_analysis_across_time``,
  ``cmh_test_across_time``, ``cmh_test``) now read the Stage-1 Parquet
  files directly.
- **Optional group filter**: ``--groups`` restricts the cache to the groups
  in this combination when the same cache file is shared across multiple
  group combinations (e.g. ``1D`` and ``AL`` caches share a timepoint).

Outputs produced (unchanged from original)
------------------------------------------
Longitudinal:
  - ``{mag}_allele_frequency_changes.tsv.gz``          — per-subject diffs
  - ``{mag}_allele_frequency_changes_no_zero-diff.tsv.gz`` — diffs excluding
    positions where the summed diff across all samples is zero (optional)
  - ``{mag}_allele_frequency_changes_mean.tsv.gz``     — mean diffs per
    (contig, gene_id, position, replicate, group)

Single:
  - ``{mag}_allele_frequency_single.tsv.gz``           — pass-through of
    per-base frequencies
  - ``{mag}_allele_frequency_no_constant.tsv.gz``      — positions where at
    least one nucleotide varies across samples (optional)

Shared constant (from ``_allele_freq_common.py``)
--------------------------------------------------
``NUCLEOTIDES`` — column names for the four per-base frequency columns
(``A_frequency``, ``T_frequency``, ``G_frequency``, ``C_frequency``).
These names are established in Stage 1 by ``calculate_frequencies()`` and
must be consistent with the diff-column names used here (``{nuc}_diff``).
"""

import argparse
import gc
import logging
import os
import sys
import time

import pandas as pd

# NUCLEOTIDES is the single source of truth for frequency column names.
# Stage 1 writes these columns; Stage 2 reads and diffs them.
from alleleflux.scripts.analysis._allele_freq_common import NUCLEOTIDES
from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cache loading  (replaces the old QC+profile reading block)
# ---------------------------------------------------------------------------

def load_cache_files(cache_paths):
    """Read one or two Stage-1 Parquet cache files into a single DataFrame.

    Previously this step was ``pd.read_csv(longitudinal_file, sep="\\t")``.
    Parquet is used here because:
    - Columnar layout allows downstream steps to load only needed columns.
    - Snappy decompression is several times faster than gzip for large files.
    - No schema inference overhead on each read (types are stored in metadata).

    Parameters
    ----------
    cache_paths : list of str
        One path for single data; two paths (one per timepoint) for
        longitudinal data.

    Returns
    -------
    pd.DataFrame
        Concatenated frame with all samples across all supplied cache files.
    """
    frames = []
    for p in cache_paths:
        logger.info(f"Reading cache file {p}")
        frames.append(pd.read_parquet(p))
    df = pd.concat(frames, ignore_index=True)
    logger.info(f"Loaded {len(df):,} rows from {len(cache_paths)} cache file(s).")
    return df


# ---------------------------------------------------------------------------
# Longitudinal helpers  (unchanged logic from original allele_freq.py)
# ---------------------------------------------------------------------------

def split_into_per_sample(df):
    """Split the concatenated cache DataFrame into per-sample frames.

    The cache file holds all samples for a timepoint; for longitudinal diff
    computation we need to process each sample independently before merging
    timepoints within a subject.

    Returns a list rather than a generator so the caller can ``del df`` and
    free memory before the list is consumed.
    """
    return [g for _, g in df.groupby("sample_id", sort=False)]


def create_data_dict(data_list):
    """Index per-sample frames by (subjectID, timepoint).

    Returns
    -------
    dict
        ``{subjectID: {timepoint: DataFrame}}``.  Used by
        ``calculate_allele_frequency_changes`` to pair the two timepoints
        for each subject.

    Raises
    ------
    ValueError
        If a single-sample frame contains more than one subjectID or
        more than one timepoint (indicates a data-integrity problem in the
        cache file).
    """
    data_dict = {}
    for df in data_list:
        subject_ids = df["subjectID"].unique()
        if len(subject_ids) != 1:
            raise ValueError(
                f"Multiple subjectIDs found in DataFrame for sample {df['sample_id'].iloc[0]}"
            )
        subjectID = subject_ids[0]
        timepoints = df["time"].unique()
        if len(timepoints) != 1:
            raise ValueError(
                f"Multiple timepoints found in DataFrame for sample {df['sample_id'].iloc[0]}"
            )
        timepoint = timepoints[0]

        if subjectID not in data_dict:
            data_dict[subjectID] = {}
        data_dict[subjectID][timepoint] = df
    return data_dict


def calculate_allele_frequency_changes(data_dict, output_dir, mag_id):
    """Compute per-position allele-frequency diffs between the two timepoints.

    For each subject that has data at both timepoints, merges the two
    timepoint frames on (subjectID, contig, gene_id, position, replicate,
    group) and subtracts timepoint-1 frequencies from timepoint-2 frequencies.

    Column naming convention for diffs (defined by ``NUCLEOTIDES``):
        ``{nuc}_diff = {nuc}_{timepoint_2} - {nuc}_{timepoint_1}``

    Also writes the intermediate ``_allele_frequency_changes.tsv.gz`` — a
    per-subject, per-position diff table (one row per subject × position).

    Parameters
    ----------
    data_dict : dict
        Output of ``create_data_dict``.
    output_dir : str
        Directory for the output file.
    mag_id : str
        Used in output filenames.

    Returns
    -------
    pd.DataFrame
        The concatenated per-subject diff table.  Passed to
        ``filter_constant_positions`` and then ``get_mean_change``.
    """
    logger.info("Identifying unique timepoints.")

    unique_timepoints = set()
    for subject_data in data_dict.values():
        unique_timepoints.update(subject_data.keys())
    if len(unique_timepoints) != 2:
        raise ValueError(
            f"Expected exactly 2 unique timepoints, found {len(unique_timepoints)}."
        )

    timepoint_1, timepoint_2 = unique_timepoints

    logger.info(
        f"Calculating change in allele frequency between {timepoint_1} and "
        f"{timepoint_2} for each position between the same subjectID."
    )

    # Identify subjects present at only one timepoint and log them as warnings.
    subjectIDs_timepoint1 = {
        s for s in data_dict if timepoint_1 in data_dict[s]
    }
    subjectIDs_timepoint2 = {
        s for s in data_dict if timepoint_2 in data_dict[s]
    }

    if diff := subjectIDs_timepoint1 - subjectIDs_timepoint2:
        logger.warning(
            f"SubjectIDs only in '{timepoint_1}' (no match at '{timepoint_2}'): {diff}"
        )
    if diff := subjectIDs_timepoint2 - subjectIDs_timepoint1:
        logger.warning(
            f"SubjectIDs only in '{timepoint_2}' (no match at '{timepoint_1}'): {diff}"
        )

    common_subjectIDs = [
        s for s in data_dict
        if timepoint_1 in data_dict[s] and timepoint_2 in data_dict[s]
    ]

    if not common_subjectIDs:
        raise ValueError(
            f"No common subjectIDs found between {timepoint_1} and {timepoint_2}."
        )
    logger.info(f"Common subjectIDs: {common_subjectIDs}")

    results = []
    for subjectID in common_subjectIDs:
        df_timepoint1 = data_dict[subjectID][timepoint_1]
        df_timepoint2 = data_dict[subjectID][timepoint_2]

        # Inner join: only positions present at both timepoints are kept.
        # Suffixes encode the timepoint in the merged column names so
        # diff calculation is unambiguous.
        merged_df = pd.merge(
            df_timepoint1,
            df_timepoint2,
            on=["subjectID", "contig", "gene_id", "position", "replicate", "group"],
            suffixes=(f"_{timepoint_1}", f"_{timepoint_2}"),
            how="inner",
        )

        if merged_df.empty:
            logger.warning(f"No matching positions found for subjectID {subjectID}.")
            continue

        # Compute frequency diff for each nucleotide.
        # NUCLEOTIDES is imported from _allele_freq_common to stay in sync
        # with the column names written by Stage 1 (calculate_frequencies).
        for nuc in NUCLEOTIDES:
            merged_df[f"{nuc}_diff"] = (
                merged_df[f"{nuc}_{timepoint_2}"] - merged_df[f"{nuc}_{timepoint_1}"]
            )

        # Combined coverage is used by downstream CMH/LMM tests.
        merged_df["total_coverage_combined"] = (
            merged_df[f"total_coverage_{timepoint_1}"]
            + merged_df[f"total_coverage_{timepoint_2}"]
        )

        columns_to_keep = (
            ["subjectID", "gene_id", "contig", "position", "replicate", "group"]
            + [f"total_coverage_{timepoint_1}", f"total_coverage_{timepoint_2}",
               "total_coverage_combined"]
            + [f"{nuc}_{timepoint_1}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_{timepoint_2}" for nuc in NUCLEOTIDES]
            + [f"{nuc}_diff" for nuc in NUCLEOTIDES]
        )
        results.append(merged_df[columns_to_keep])

    if not results:
        logger.error("No allele frequency changes calculated.")
        sys.exit(42)

    allele_changes = pd.concat(results, ignore_index=True)

    allele_changes.to_csv(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_changes.tsv.gz"),
        sep="\t",
        index=False,
        compression="gzip",
    )
    logger.info(
        f"Allele frequency changes saved to "
        f"{output_dir}/{mag_id}_allele_frequency_changes.tsv.gz"
    )
    return allele_changes


def filter_constant_positions(allele_df, output_dir, mag_id, data_type):
    """Remove uninformative positions before mean aggregation.

    For single data
    ---------------
    Drops positions where the frequency of every nucleotide is identical
    across all samples (no variation to test statistically).

    For longitudinal data
    ---------------------
    Drops "zero-diff" positions: those where the sum of per-subject
    frequency diffs is zero for ALL four nucleotides across all samples.
    These positions show no net change in the population and add noise to
    statistical tests.

    Parameters
    ----------
    allele_df : pd.DataFrame
        For longitudinal: the output of ``calculate_allele_frequency_changes``.
        For single: the per-base frequency frame from the cache.
    output_dir : str
        Output directory.
    mag_id : str
        Used in output filenames.
    data_type : str
        ``"single"`` or ``"longitudinal"``.

    Returns
    -------
    pd.DataFrame or None
        For longitudinal: the filtered diff DataFrame (passed to
        ``get_mean_change``).
        For single: ``None`` (output written to disk inside this function).
    """
    grouped_df = allele_df.groupby(["contig", "position"], dropna=False)

    if data_type == "single":
        logger.info(
            "Identifying positions where allele frequency values across all "
            "samples for each nucleotide is constant"
        )
        # nunique > 1 for any nucleotide means the position is variable.
        agg = grouped_df[NUCLEOTIDES].nunique()
        groups_to_keep = (agg > 1).any(axis=1)

        positions_kept = groups_to_keep.sum()
        total_positions = agg.shape[0]
        positions_removed = total_positions - positions_kept
        logger.info(f"Found {positions_removed:,} constant positions. Filtering...")

        keep_index = groups_to_keep[groups_to_keep].index
        filtered_df = (
            allele_df.set_index(["contig", "position"]).loc[keep_index].reset_index()
        )
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, "
            f"Positions removed: {positions_removed:,}"
        )

        filtered_df.to_csv(
            os.path.join(output_dir, f"{mag_id}_allele_frequency_no_constant.tsv.gz"),
            sep="\t",
            index=False,
            compression="gzip",
        )
        logger.info(
            f"Allele frequency changes with no constant positions saved to "
            f"{output_dir}/{mag_id}_allele_frequency_no_constant.tsv.gz"
        )
        return None

    elif data_type == "longitudinal":
        logger.info(
            "Identifying positions where sum of the difference values across "
            "all samples for all nucleotides is zero, called zero-diff positions"
        )

        # Build diff column names from NUCLEOTIDES (e.g. "A_frequency_diff").
        diff_cols = [f"{nuc}_diff" for nuc in NUCLEOTIDES]

        grouped_sums = grouped_df[diff_cols].sum()

        # A position is zero-diff only if ALL four sums are exactly zero.
        is_all_zero = (
            (grouped_sums["A_frequency_diff"] == 0)
            & (grouped_sums["T_frequency_diff"] == 0)
            & (grouped_sums["G_frequency_diff"] == 0)
            & (grouped_sums["C_frequency_diff"] == 0)
        )
        zero_positions = is_all_zero[is_all_zero].index
        logger.info(f"Found {len(zero_positions):,} zero-diff positions.")

        logger.info("Filtering zero-diff positions")
        ac_indexed = allele_df.set_index(["contig", "position"], drop=False)
        keep_mask_ac = ~ac_indexed.index.isin(zero_positions)
        filtered_allele_changes = ac_indexed[keep_mask_ac].copy()
        filtered_allele_changes.reset_index(drop=True, inplace=True)

        total_positions = grouped_sums.shape[0]
        positions_removed = len(zero_positions)
        positions_kept = total_positions - positions_removed
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, "
            f"Positions removed: {positions_removed:,}"
        )

        filtered_allele_changes.to_csv(
            os.path.join(
                output_dir, f"{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
            ),
            sep="\t",
            index=False,
            compression="gzip",
        )
        logger.info(
            f"Allele frequency changes with no zero diff positions saved to "
            f"{output_dir}/{mag_id}_allele_frequency_changes_no_zero-diff.tsv.gz"
        )
        return filtered_allele_changes


def get_mean_change(allele_changes, mag_id, output_dir):
    """Aggregate per-subject diffs into per-(position, group, replicate) means.

    Input rows are per-subject (one row per subjectID × position × group ×
    replicate).  This function collapses the subjectID dimension by computing
    the mean diff for each nucleotide, producing the ``_changes_mean.tsv.gz``
    file consumed by single-sample tests, LMM, CMH between-groups, and
    preprocessing rules.

    The ``subjectID`` column is replaced by ``subjectID_count`` (nunique) so
    downstream code knows how many subjects contributed to each mean.

    Parameters
    ----------
    allele_changes : pd.DataFrame
        Optionally zero-diff-filtered output of
        ``calculate_allele_frequency_changes`` /
        ``filter_constant_positions``.
    mag_id : str
        Used in output filename.
    output_dir : str
        Output directory.

    Returns
    -------
    pd.DataFrame
        Mean-aggregated frame (also written to disk).
    """
    logger.info(
        "Calculating mean changes in allele frequencies for subjectIDs "
        "present in the same replicate and group."
    )
    # Mean of each diff column; count of unique subjects per group.
    agg_dict = {f"{nuc}_diff": "mean" for nuc in NUCLEOTIDES}
    agg_dict["subjectID"] = "nunique"

    mean_changes_df = (
        allele_changes.groupby(
            ["contig", "gene_id", "position", "replicate", "group"], dropna=False
        )
        .agg(agg_dict)
        .reset_index()
    )

    mean_changes_df.rename(
        columns={f"{nuc}_diff": f"{nuc}_diff_mean" for nuc in NUCLEOTIDES},
        inplace=True,
    )
    mean_changes_df.rename(columns={"subjectID": "subjectID_count"}, inplace=True)

    mean_changes_df.to_csv(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_changes_mean.tsv.gz"),
        index=False,
        sep="\t",
        compression="gzip",
    )
    logger.info(
        f"Mean change in allele frequency for MAG {mag_id} saved to {output_dir}"
    )
    return mean_changes_df


# ---------------------------------------------------------------------------
# Single-data path
# ---------------------------------------------------------------------------

def process_single_data(allele_df, output_dir, mag_id, disable_filtering):
    """Write the single-data outputs (pass-through + optional constant filter)."""
    output_fpath = os.path.join(output_dir, f"{mag_id}_allele_frequency_single.tsv.gz")
    logger.info(f"Writing allele frequencies (single data) to {output_fpath}")
    allele_df.to_csv(output_fpath, sep="\t", index=False, compression="gzip")

    if not disable_filtering:
        logger.info("Filtering constant allele frequency positions (single data).")
        filter_constant_positions(allele_df, output_dir, mag_id, data_type="single")
    else:
        logger.info("User disabled filtering of constant positions.")


# ---------------------------------------------------------------------------
# Longitudinal data path
# ---------------------------------------------------------------------------

def process_longitudinal_data(allele_df, output_dir, mag_id, disable_filtering):
    """Orchestrate the longitudinal diff pipeline.

    Steps
    -----
    1. Split the concatenated cache frame back into per-sample frames
       (``split_into_per_sample``).
    2. Index frames by (subjectID, timepoint) for pairwise diff
       (``create_data_dict``).
    3. Compute per-position frequency diffs between the two timepoints
       (``calculate_allele_frequency_changes``).
    4. Optionally filter zero-diff positions (``filter_constant_positions``).
    5. Compute per-(position, group, replicate) mean diffs
       (``get_mean_change``).

    Memory management: intermediate DataFrames are deleted and
    ``gc.collect()`` is called at each step to keep the peak footprint low.
    The cache holds both timepoints in one frame, which is larger than the
    old single-timepoint longitudinal file, so explicit cleanup is important.
    """
    data_list = split_into_per_sample(allele_df)
    del allele_df   # free the large concatenated cache frame
    gc.collect()

    data_dict = create_data_dict(data_list)
    del data_list
    gc.collect()

    allele_changes = calculate_allele_frequency_changes(data_dict, output_dir, mag_id)

    if not disable_filtering:
        logger.info("Filtering zero-diff positions for longitudinal data.")
        allele_changes = filter_constant_positions(
            allele_changes, output_dir, mag_id, data_type="longitudinal"
        )
    else:
        logger.info("User disabled zero-diff filtering for longitudinal data.")

    get_mean_change(allele_changes, mag_id, output_dir)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description=(
            "Stage 2: Compute allele-frequency diffs for a MAG from one or two "
            "per-timepoint Parquet cache files produced by allele_freq_cache.py "
            "(alleleflux-cache-allele-freq)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--magID",
        required=True,
        type=str,
        help="MAG ID to process.",
    )
    parser.add_argument(
        "--cache_files",
        required=True,
        nargs="+",
        type=str,
        help=(
            "Path(s) to per-(MAG, timepoint) Parquet cache file(s) written by "
            "allele_freq_cache.py.  Provide two paths for longitudinal data "
            "(one per timepoint), one path for single data."
        ),
    )
    parser.add_argument(
        "--data_type",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        type=str,
        default=None,
        help=(
            "Optional list of group names to retain before computing diffs.  "
            "Use when the cache file covers multiple group combinations sharing "
            "a timepoint; pass only the groups relevant to this combination."
        ),
    )
    parser.add_argument(
        "--disable_zero_diff_filtering",
        dest="disable_filtering",
        action="store_true",
        default=False,
        help=(
            "Skip the constant/zero-diff filtering step.  "
            "For single: skip removing positions where all nucleotide frequencies "
            "are constant.  For longitudinal: skip removing positions where the "
            "sum of frequency differences across all samples is zero."
        ),
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="Path to output directory.",
    )

    args = parser.parse_args()
    start_time = time.time()
    mag_id = args.magID

    # Validate cache file count against data type.
    if args.data_type == "longitudinal" and len(args.cache_files) != 2:
        parser.error(
            f"Longitudinal data_type requires exactly 2 cache files, "
            f"got {len(args.cache_files)}."
        )
    if args.data_type == "single" and len(args.cache_files) != 1:
        parser.error(
            f"Single data_type requires exactly 1 cache file, "
            f"got {len(args.cache_files)}."
        )

    # ------------------------------------------------------------------
    # 1. Load cache files
    # ------------------------------------------------------------------
    # For longitudinal data: the two files hold different timepoints.
    # Both are concatenated; ``create_data_dict`` later splits them back
    # by (subjectID, time) for pairwise diff computation.
    allele_df = load_cache_files(args.cache_files)

    # ------------------------------------------------------------------
    # 2. Optionally restrict to the groups in this combination
    # ------------------------------------------------------------------
    # A Parquet cache file may hold samples from multiple groups (because
    # different group combinations share the same timepoint).  When this
    # combination only involves a subset of groups, filter before diffing
    # to avoid cross-group contamination.
    if args.groups:
        before = len(allele_df)
        allele_df = allele_df[allele_df["group"].isin(args.groups)].copy()
        logger.info(
            f"Filtered cache to groups {args.groups}: "
            f"{before:,} -> {len(allele_df):,} rows."
        )
        if allele_df.empty:
            logger.error(
                f"No rows remain after filtering to groups {args.groups}. Exiting."
            )
            sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 3. Branch on data type
    # ------------------------------------------------------------------
    if args.data_type == "single":
        process_single_data(allele_df, args.output_dir, mag_id, args.disable_filtering)
    else:
        process_longitudinal_data(
            allele_df, args.output_dir, mag_id, args.disable_filtering
        )

    logger.info(f"Total time taken: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
