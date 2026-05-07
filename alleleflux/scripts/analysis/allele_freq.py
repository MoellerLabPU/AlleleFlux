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

Outputs produced
----------------
Longitudinal:
  - ``{mag}_allele_frequency_changes.parquet``         — per-subject diffs
    (Parquet/Snappy; archival only — not consumed by downstream rules)
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
from concurrent.futures import ThreadPoolExecutor

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

    When two paths are supplied (longitudinal), both files are read
    concurrently via ThreadPoolExecutor — Parquet I/O releases the GIL so
    the reads overlap in wall-clock time.

    Parameters
    ----------
    cache_paths : list of str
        One path for single data; two paths (one per timepoint) for
        longitudinal data.

    Returns
    -------
    pd.DataFrame
        Concatenated frame with all samples across all supplied cache files.
        For longitudinal data the first cache file's timepoint appears first
        in the frame, preserving the tp1/tp2 ordering used by
        ``calculate_allele_frequency_changes``.
    """
    def _read(p):
        logger.info(f"Reading cache file {p}")
        return pd.read_parquet(p)

    if len(cache_paths) == 1:
        frames = [_read(cache_paths[0])]
    else:
        with ThreadPoolExecutor(max_workers=len(cache_paths)) as ex:
            frames = list(ex.map(_read, cache_paths))

    df = pd.concat(frames, ignore_index=True)
    logger.info(f"Loaded {len(df):,} rows from {len(cache_paths)} cache file(s).")
    return df


# ---------------------------------------------------------------------------
# Longitudinal helpers
# ---------------------------------------------------------------------------

def calculate_allele_frequency_changes(allele_df, output_dir, mag_id):
    """Compute per-position allele-frequency diffs via a single vectorized merge.

    Replaces the original per-subject loop (N sequential ``pd.merge`` calls)
    with one global inner join on
    (subjectID, contig, gene_id, position, replicate, group).  The subjectID
    key prevents cross-subject pairing, so the result is identical to the
    loop but computed in a single C-level hash join.

    Timepoint ordering is preserved from the cache-file concatenation order:
    the first cache file supplied to ``load_cache_files`` becomes
    ``timepoint_1`` (reference), the second becomes ``timepoint_2`` (later).
    Snakemake passes cache files in tp1/tp2 order, so this matches the
    combination label (e.g. ``5mo_10mo``).

    Column naming convention for diffs (defined by ``NUCLEOTIDES``):
        ``{nuc}_diff = {nuc}_{timepoint_2} - {nuc}_{timepoint_1}``

    Writes ``{mag}_allele_frequency_changes.parquet`` (Parquet/Snappy) as an
    archival artifact.  No downstream Snakemake rule reads this file; Parquet
    is used because it is significantly faster to write than gzip TSV for a
    table of this width and row count.

    Parameters
    ----------
    allele_df : pd.DataFrame
        Concatenated output of ``load_cache_files`` (both timepoints).
    output_dir : str
        Directory for the output file.
    mag_id : str
        Used in output filenames.

    Returns
    -------
    pd.DataFrame
        Per-subject, per-position diff table passed to
        ``filter_constant_positions`` and then ``get_mean_change``.
    """
    # unique() preserves first-occurrence order, so timepoint_1 is the
    # timepoint from the first cache file (the earlier / reference timepoint).
    timepoints = allele_df["time"].unique()
    if len(timepoints) != 2:
        raise ValueError(
            f"Expected exactly 2 unique timepoints, found {len(timepoints)}."
        )
    timepoint_1, timepoint_2 = timepoints
    logger.info(
        f"Calculating allele frequency changes between {timepoint_1} and {timepoint_2}."
    )

    drop_cols = ["time", "sample_id"]
    df_tp1 = allele_df[allele_df["time"] == timepoint_1].drop(columns=drop_cols, errors="raise")
    df_tp2 = allele_df[allele_df["time"] == timepoint_2].drop(columns=drop_cols, errors="raise")

    subjects_tp1 = set(df_tp1["subjectID"].unique())
    subjects_tp2 = set(df_tp2["subjectID"].unique())
    if diff := subjects_tp1 - subjects_tp2:
        logger.warning(
            f"SubjectIDs only in '{timepoint_1}' (no match at '{timepoint_2}'): {diff}"
        )
    if diff := subjects_tp2 - subjects_tp1:
        logger.warning(
            f"SubjectIDs only in '{timepoint_2}' (no match at '{timepoint_1}'): {diff}"
        )

    merged_df = pd.merge(
        df_tp1,
        df_tp2,
        on=["subjectID", "contig", "gene_id", "position", "replicate", "group"],
        suffixes=(f"_{timepoint_1}", f"_{timepoint_2}"),
        how="inner",
    )
    del df_tp1, df_tp2
    gc.collect()

    if merged_df.empty:
        logger.error(
            f"No matching positions found across timepoints {timepoint_1} and {timepoint_2}."
        )
        sys.exit(42)

    for nuc in NUCLEOTIDES:
        merged_df[f"{nuc}_diff"] = (
            merged_df[f"{nuc}_{timepoint_2}"] - merged_df[f"{nuc}_{timepoint_1}"]
        )
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
    allele_changes = merged_df[columns_to_keep].copy()
    del merged_df
    gc.collect()

    allele_changes.to_parquet(
        os.path.join(output_dir, f"{mag_id}_allele_frequency_changes.parquet"),
        engine="pyarrow",
        compression="snappy",
        index=False,
    )
    logger.info(
        f"Allele frequency changes saved to "
        f"{output_dir}/{mag_id}_allele_frequency_changes.parquet"
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
    1. Compute per-position frequency diffs via a single vectorized merge
       (``calculate_allele_frequency_changes``).
    2. Optionally filter zero-diff positions (``filter_constant_positions``).
    3. Compute per-(position, group, replicate) mean diffs
       (``get_mean_change``).
    """
    allele_changes = calculate_allele_frequency_changes(allele_df, output_dir, mag_id)
    del allele_df
    gc.collect()

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
    # Both are concatenated; ``calculate_allele_frequency_changes`` later
    # splits by timepoint and performs a single vectorized merge.
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
        allele_df["group"] = allele_df["group"].astype(str)
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
