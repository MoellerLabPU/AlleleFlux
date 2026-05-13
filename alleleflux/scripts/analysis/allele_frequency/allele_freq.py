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
  - ``{mag}_allele_frequency_changes.parquet``             — per-subject diffs
    (Parquet/Snappy; archival only — not consumed by downstream rules)
  - ``{mag}_allele_frequency_changes_no_zero-diff.parquet`` — diffs excluding
    positions where the summed diff across all samples is zero (optional)
  - ``{mag}_allele_frequency_changes_mean.parquet``        — mean diffs per
    (contig, gene_id, position, replicate, group)

Single:
  - ``{mag}_allele_frequency_single.parquet``          — pass-through of
    per-base frequencies
  - ``{mag}_allele_frequency_no_constant.parquet``     — positions where at
    least one nucleotide varies across samples (optional)

All outputs use Parquet/Snappy compression. Downstream consumers read them
via ``load_allele_freq_inputs()`` in utilities.py which handles both .parquet
and .tsv.gz by extension — no consumer-side changes needed.

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
from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import NUCLEOTIDES
from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cache loading  (replaces the old QC+profile reading block)
# ---------------------------------------------------------------------------

# Columns Stage 2 reads from the cache.
#
# Longitudinal: only the frequency columns are needed.  Raw A/C/G/T counts
# and sample_id are skipped via ``columns=`` (pyarrow pushes the selection
# down to the parquet reader, so unread columns never enter pandas).  This
# is roughly a 30% reduction in both I/O and peak RAM vs reading every
# column.  CMH/LMM bypass Stage 2 and read the cache files directly when
# they need raw counts, so dropping them here does not affect those rules.
#
# Single: keep raw A/C/G/T because single-data CMH and other downstream
# consumers read the Stage 2 output ({mag}_allele_frequency_single.parquet
# or {mag}_allele_frequency_no_constant.parquet) and expect raw counts to
# be present.  sample_id is still dropped — it's not consumed anywhere
# downstream and is one of the heaviest string columns in the cache.
_STAGE2_COLUMNS_LONGITUDINAL = [
    "contig", "position", "gene_id",
    "total_coverage",
    "group", "subjectID", "replicate",
    "A_frequency", "T_frequency", "G_frequency", "C_frequency",
    "time",
]
_STAGE2_COLUMNS_SINGLE = [
    "contig", "position", "gene_id",
    "total_coverage",
    "A", "C", "G", "T",
    "group", "subjectID", "replicate",
    "A_frequency", "T_frequency", "G_frequency", "C_frequency",
]


def load_cache_files(cache_paths, data_type):
    """Read one or two Stage-1 Parquet cache files into a single DataFrame.

    Only the columns Stage 2 needs are read from disk via ``columns=``.
    pyarrow pushes the column selection down to the parquet reader so
    unread columns never enter memory — typically a 25–30 % reduction in
    both I/O and peak pandas footprint.

    When two paths are supplied (longitudinal), both files are read
    concurrently via ThreadPoolExecutor — Parquet I/O releases the GIL so
    the reads overlap in wall-clock time.

    Parameters
    ----------
    cache_paths : list of str
        One path for single data; two paths (one per timepoint) for
        longitudinal data.
    data_type : str
        ``"single"`` or ``"longitudinal"``.  Selects the column keep-list:
        longitudinal drops A/C/G/T (raw counts unused by diff computation),
        single keeps them (downstream CMH expects them in the output).

    Returns
    -------
    pd.DataFrame
        Concatenated frame with all samples across all supplied cache files.
        For longitudinal data the first cache file's timepoint appears first
        in the frame, preserving the tp1/tp2 ordering used by
        ``calculate_allele_frequency_changes``.
    """
    columns = (
        _STAGE2_COLUMNS_LONGITUDINAL
        if data_type == "longitudinal"
        else _STAGE2_COLUMNS_SINGLE
    )

    def _read(p):
        logger.info(f"Reading cache file {p}")
        return pd.read_parquet(p, columns=columns)

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

def calculate_allele_frequency_changes(allele_df, output_dir, mag_id, save_archival=False):
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

    Memory shape
    ------------
    Peak memory inside this function is the merged frame (tp1 + tp2 columns
    side-by-side + diff columns).  We always drop the three
    ``total_coverage_*`` columns after the merge — they are only consumed by
    the archival output.  The per-timepoint frequency columns
    (``A_frequency_{tp1}``, etc.) are kept in the returned frame because
    ``lmm_analysis_across_time`` reads them out of the
    ``_changes_no_zero-diff.parquet`` file in "suffix" mode.

    Parameters
    ----------
    allele_df : pd.DataFrame
        Concatenated output of ``load_cache_files`` (both timepoints).
    output_dir : str
        Directory for the output file (only used when ``save_archival``).
    mag_id : str
        Used in output filenames.
    save_archival : bool, optional
        When True, write ``{mag}_allele_frequency_changes.parquet`` — a
        per-subject diff table that includes the original tp1/tp2
        frequency and coverage columns.  No downstream Snakemake rule
        reads this file, so default is False to save memory and disk.

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

    # ``sample_id`` is no longer loaded from the cache (see load_cache_files),
    # so only ``time`` needs to be stripped before the merge.
    df_tp1 = allele_df[allele_df["time"] == timepoint_1].drop(columns=["time"])
    df_tp2 = allele_df[allele_df["time"] == timepoint_2].drop(columns=["time"])

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

    # Compute diffs.  These are the only post-merge columns consumed by
    # downstream steps (filter_constant_positions, get_mean_change).
    for nuc in NUCLEOTIDES:
        merged_df[f"{nuc}_diff"] = (
            merged_df[f"{nuc}_{timepoint_2}"] - merged_df[f"{nuc}_{timepoint_1}"]
        )

    diff_cols = [f"{nuc}_diff" for nuc in NUCLEOTIDES]
    key_cols = ["subjectID", "gene_id", "contig", "position", "replicate", "group"]
    tp1_freq_cols = [f"{nuc}_{timepoint_1}" for nuc in NUCLEOTIDES]
    tp2_freq_cols = [f"{nuc}_{timepoint_2}" for nuc in NUCLEOTIDES]

    if save_archival:
        # Archival path — write the full per-subject diff table including
        # total_coverage columns.  Not consumed by any downstream Snakemake
        # rule, hence opt-in.  Computing total_coverage_combined here keeps
        # the column out of the post-merge memory peak when archival is off.
        merged_df["total_coverage_combined"] = (
            merged_df[f"total_coverage_{timepoint_1}"]
            + merged_df[f"total_coverage_{timepoint_2}"]
        )
        archival_cols = (
            key_cols
            + [f"total_coverage_{timepoint_1}", f"total_coverage_{timepoint_2}",
               "total_coverage_combined"]
            + tp1_freq_cols + tp2_freq_cols + diff_cols
        )
        merged_df[archival_cols].to_parquet(
            os.path.join(output_dir, f"{mag_id}_allele_frequency_changes.parquet"),
            engine="pyarrow",
            compression="snappy",
            index=False,
        )
        logger.info(
            f"Allele frequency changes saved to "
            f"{output_dir}/{mag_id}_allele_frequency_changes.parquet"
        )

    # Keep tp1/tp2 frequency columns for LMM across-time (which reads them
    # out of {mag}_allele_frequency_changes_no_zero-diff.parquet in "suffix"
    # mode), plus the diffs that get_mean_change aggregates.  Drop the three
    # total_coverage_* columns — they're archival-only and freeing them
    # here saves ~24 B/row before the next groupby.
    keep_cols = key_cols + tp1_freq_cols + tp2_freq_cols + diff_cols
    allele_changes = merged_df[keep_cols]
    del merged_df
    gc.collect()
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
    grouped_df = allele_df.groupby(["contig", "position"], dropna=False, observed=True)

    # Build a (contig, position) MultiIndex directly from the two key columns
    # instead of going through ``allele_df.set_index([...])``.  ``set_index``
    # rebuilds the whole DataFrame (every column block gets a new owner) just
    # to expose a MultiIndex we throw away after ``.isin()``.  Constructing
    # the MultiIndex from numpy views of the key columns avoids touching the
    # rest of the frame.
    def _position_index(df):
        return pd.MultiIndex.from_arrays(
            [df["contig"].to_numpy(), df["position"].to_numpy()],
            names=["contig", "position"],
        )

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

        keep_positions = groups_to_keep[groups_to_keep].index
        mask = _position_index(allele_df).isin(keep_positions)
        filtered_df = allele_df[mask].reset_index(drop=True)
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, "
            f"Positions removed: {positions_removed:,}"
        )

        out_path = os.path.join(output_dir, f"{mag_id}_allele_frequency_no_constant.parquet")
        filtered_df.to_parquet(out_path, compression="snappy", index=False)
        logger.info(
            f"Allele frequency (no constant positions) saved to {out_path}"
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
        keep_mask_ac = ~_position_index(allele_df).isin(zero_positions)
        filtered_allele_changes = allele_df[keep_mask_ac].reset_index(drop=True)

        total_positions = grouped_sums.shape[0]
        positions_removed = len(zero_positions)
        positions_kept = total_positions - positions_removed
        logger.info(
            f"Total positions: {total_positions:,}, Positions kept: {positions_kept:,}, "
            f"Positions removed: {positions_removed:,}"
        )

        out_path = os.path.join(
            output_dir, f"{mag_id}_allele_frequency_changes_no_zero-diff.parquet"
        )
        filtered_allele_changes.to_parquet(out_path, compression="snappy", index=False)
        logger.info(
            f"Allele frequency changes (no zero-diff) saved to {out_path}"
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

    # observed=True avoids the PerformanceWarning on categorical groupby keys
    # and is significantly faster when group/replicate are categorical.
    mean_changes_df = (
        allele_changes.groupby(
            ["contig", "gene_id", "position", "replicate", "group"],
            dropna=False,
            observed=True,
        )
        .agg(agg_dict)
        .reset_index()
    )

    mean_changes_df.rename(
        columns={f"{nuc}_diff": f"{nuc}_diff_mean" for nuc in NUCLEOTIDES},
        inplace=True,
    )
    mean_changes_df.rename(columns={"subjectID": "subjectID_count"}, inplace=True)

    out_path = os.path.join(output_dir, f"{mag_id}_allele_frequency_changes_mean.parquet")
    mean_changes_df.to_parquet(out_path, compression="snappy", index=False)
    logger.info(
        f"Mean allele frequency changes for MAG {mag_id} saved to {out_path}"
    )
    return mean_changes_df


# ---------------------------------------------------------------------------
# Single-data path
# ---------------------------------------------------------------------------

def process_single_data(allele_df, output_dir, mag_id, disable_filtering):
    """Write the single-data outputs (pass-through + optional constant filter)."""
    output_fpath = os.path.join(output_dir, f"{mag_id}_allele_frequency_single.parquet")
    logger.info(f"Writing allele frequencies (single data) to {output_fpath}")
    allele_df.to_parquet(output_fpath, compression="snappy", index=False)

    if not disable_filtering:
        logger.info("Filtering constant allele frequency positions (single data).")
        filter_constant_positions(allele_df, output_dir, mag_id, data_type="single")
    else:
        logger.info("User disabled filtering of constant positions.")


# ---------------------------------------------------------------------------
# Longitudinal data path
# ---------------------------------------------------------------------------

def process_longitudinal_data(allele_df, output_dir, mag_id, disable_filtering, save_archival=False):
    """Orchestrate the longitudinal diff pipeline.

    Steps
    -----
    1. Compute per-position frequency diffs via a single vectorized merge
       (``calculate_allele_frequency_changes``).
    2. Optionally filter zero-diff positions (``filter_constant_positions``).
    3. Compute per-(position, group, replicate) mean diffs
       (``get_mean_change``).

    ``save_archival`` forwards through to ``calculate_allele_frequency_changes``;
    when False (the default) the per-subject diff parquet is not written and
    its column slab is freed before the next groupby.
    """
    allele_changes = calculate_allele_frequency_changes(
        allele_df, output_dir, mag_id, save_archival=save_archival
    )
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
        "--save_archival_changes",
        action="store_true",
        default=False,
        help=(
            "Write {mag}_allele_frequency_changes.parquet — the per-subject diff "
            "table including the original tp1/tp2 frequency and coverage columns.  "
            "Not consumed by any downstream Snakemake rule; off by default to "
            "save ~150 GB of transient memory and one large disk write per MAG."
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
    # ``data_type`` controls which columns are loaded (longitudinal drops
    # raw A/C/G/T counts; see _STAGE2_COLUMNS_* above).
    allele_df = load_cache_files(args.cache_files, args.data_type)

    # ------------------------------------------------------------------
    # 2. Optionally restrict to the groups in this combination
    # ------------------------------------------------------------------
    # A Parquet cache file may hold samples from multiple groups (because
    # different group combinations share the same timepoint).  When this
    # combination only involves a subset of groups, filter before diffing
    # to avoid cross-group contamination.
    if args.groups:
        before = len(allele_df)
        # Compare as strings so a categorical group column (from the new
        # cache format) and numeric group values both match cleanly against
        # args.groups (which is a list of strings from the CLI).  We cast a
        # *view* via .astype(str) — the original column remains untouched
        # so we can reapply the categorical dtype below without re-inferring
        # from a fresh string column.
        keep_mask = allele_df["group"].astype(str).isin(args.groups)
        allele_df = allele_df[keep_mask].reset_index(drop=True)
        logger.info(
            f"Filtered cache to groups {args.groups}: "
            f"{before:,} -> {len(allele_df):,} rows."
        )
        if allele_df.empty:
            logger.error(
                f"No rows remain after filtering to groups {args.groups}. Exiting."
            )
            sys.exit(1)

    # Apply / refresh categorical dtypes.  The new cache format already stores
    # these as dictionary-encoded (read back as categorical), but we still:
    #   - convert any object columns that arrived from older cache files; and
    #   - call remove_unused_categories() so the group filter above doesn't
    #     leave dangling levels that inflate every downstream groupby.
    # contig/gene_id are deliberately excluded — too many unique values for
    # categorical to pay off in pandas, even though the parquet file
    # dictionary-encodes them on disk.
    for col in ("group", "time", "subjectID", "replicate"):
        if col in allele_df.columns:
            s = allele_df[col]
            if not isinstance(s.dtype, pd.CategoricalDtype):
                allele_df[col] = s.astype("category")
            else:
                allele_df[col] = s.cat.remove_unused_categories()

    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # 3. Branch on data type
    # ------------------------------------------------------------------
    if args.data_type == "single":
        process_single_data(allele_df, args.output_dir, mag_id, args.disable_filtering)
    else:
        process_longitudinal_data(
            allele_df,
            args.output_dir,
            mag_id,
            args.disable_filtering,
            save_archival=args.save_archival_changes,
        )

    logger.info(f"Total time taken: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
