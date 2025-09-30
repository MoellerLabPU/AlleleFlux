#!/usr/bin/env python3
"""
Compute mean per-position coverage across samples for MAGs.

This script processes all *_metadata.tsv files in a directory, where each metadata file
lists sample_id and file_path to each sample's profiled MAG file. It aggregates total_coverage
for each (contig, position) across all samples and reports mean coverage per position.
By default, averages are calculated only over samples where a position is present (absent positions are not treated as zero).
Use --include_zeros to treat absent positions as zero when computing means and stds.

Usage:
    alleleflux-mean-coverage --rootDir /path/to/metadata_files --output_dir /path/to/output

Output:
A TSV per MAG with columns:
- contig, position, mean_coverage, n_present, n_samples
- mean_coverage_group_<group> (if group column present)
- mean_coverage_group_<group>_time_<time> (if both group and time present)
- n_present_group_<group> and n_present_group_<group>_time_<time>
    (number of samples used for that mean; respects include_zeros mode)
- std_coverage (overall population std, ddof=0)
- std_coverage_group_<group> and std_coverage_group_<group>_time_<time>

Additionally, this script reports allele-frequency statistics (ATGC only) from the per-sample allele counts:
- mean_freq_<ALLELE> and std_freq_<ALLELE> (overall), where freq = <ALLELE> / total_coverage per sample and <ALLELE> in {A,C,G,T}
- Grouped counterparts for each present group/time combination, e.g., mean_freq_T_group_<g> and mean_freq_T_group_<g>_time_<t>

Important: Coverage means/stds honor the include_zeros option (i.e., absent positions can be treated as 0 coverage when averaging). However, allele frequencies are always averaged only over samples with positive coverage at that position (present-only), regardless of include_zeros, to avoid biasing frequencies toward zero.

Notes:
- Expects profile files to have columns: contig, position, total_coverage and allele columns A,C,G,T
- Duplicate (contig, position) within a file will raise an error
- Output sorted by contig, position
"""

import argparse
import logging
import multiprocessing
import os
import re
from collections import defaultdict
from glob import glob
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# Constants
ALLELES = ["A", "C", "G", "T"]


def read_mag_metadata(metadata_path: Path) -> pd.DataFrame:
    """
    Read and validate MAG metadata file containing sample information.

    The metadata file must be a TSV with at least 'sample_id' and 'file_path' columns.
    Optional columns include 'group' and 'time' for grouping samples.

    Args:
        metadata_path: Path to the TSV metadata file.

    Returns:
        A pandas DataFrame with normalized string columns.

    Raises:
        ValueError: If required columns are missing or the file is empty.
    """
    # Load the metadata from TSV
    df = pd.read_csv(metadata_path, sep="\t")
    # Check for required columns
    required = {"sample_id", "file_path"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in MAG metadata '{metadata_path}': {sorted(missing)}"
        )
    # Normalize string types by stripping whitespace
    df["sample_id"] = df["sample_id"].astype(str).str.strip()
    df["file_path"] = df["file_path"].astype(str).str.strip()
    # Normalize optional grouping columns if present
    if "group" in df.columns:
        df["group"] = df["group"].astype(str).str.strip()
    if "time" in df.columns:
        df["time"] = df["time"].astype(str).str.strip()
    # Ensure not empty
    if df.empty:
        raise ValueError(f"No rows in MAG metadata file: {metadata_path}")
    return df


def read_profile(sample_file: str) -> pd.DataFrame:
    """
    Read one sample's coverage profile and return per-position coverage data.

    The profile file should be a TSV with columns: contig, position, total_coverage
    and allele count columns A,C,G,T for downstream allele-frequency statistics.

    Args:
        sample_file: Path to the sample's profile TSV file.

    Returns:
        A DataFrame with columns ['contig', 'position', 'total_coverage', 'A', 'C', 'G', 'T'].
        If the file has no data, returns an empty DataFrame with the expected columns.

    Raises:
        ValueError: If duplicate (contig, position) pairs exist in the file.
    """
    base_cols = ["contig", "position", "total_coverage"]
    usecols = base_cols + ALLELES

    logger.debug(f"Reading sample profile: {sample_file}")
    dtype = {"contig": str, "position": int, "total_coverage": float}
    dtype.update({allele: float for allele in ALLELES})

    df = pd.read_csv(sample_file, sep="\t", usecols=usecols, dtype=dtype)

    if df.empty:
        logger.warning(
            f"Profile file not found; sample excluded (counts as zero only if --include_zeros): {sample_file}"
        )
        return pd.DataFrame(columns=usecols)

    dup = df.duplicated(subset=["contig", "position"], keep=False)
    if dup.any():
        raise ValueError(
            f"Duplicate (contig, position) pairs found in {sample_file}:\n"
            f"{df[dup].sort_values(['contig', 'position'])}"
        )
    return df


def load_and_combine_sample_data(meta: pd.DataFrame) -> pd.DataFrame:
    """
    Load and combine coverage data from all sample profile files listed in metadata.

    Args:
        meta: DataFrame with sample metadata including file_path column

    Returns:
        Combined DataFrame with all sample data tagged with sample_id, group, time
    """
    all_data = []
    # Iterate through each sample in the metadata
    for _, row in meta.iterrows():
        # Extract file path for this sample's profile
        fpath = Path(row["file_path"])
        # Skip and warn if the profile file does not exist (treat as all zeros)
        if not fpath.exists():
            logger.warning(f"Profile file not found; treated as all zeros: {fpath}")
            continue
        # Read per-position coverage data for this sample
        df = read_profile(fpath)
        # Skip if no coverage data in this sample
        if df.empty:
            continue
        # Tag the data with sample identifier
        df["sample_id"] = row["sample_id"]
        # Tag with group info if available in metadata
        if "group" in meta.columns:
            df["group"] = row["group"]
        # Tag with time info if available
        if "time" in meta.columns:
            df["time"] = row["time"]
        # Add this sample's data to the collection
        all_data.append(df)

    if not all_data:
        logger.warning("No coverage data found in any sample.")
        return pd.DataFrame()

    # Concatenate all sample DataFrames into a single comprehensive DataFrame
    return pd.concat(all_data, ignore_index=True)


def compute_allele_frequencies(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute allele frequencies and their squares for statistical calculations.

    Args:
        df: DataFrame with allele count columns A, C, G, T

    Returns:
        Modified DataFrame with frequency columns added
    """
    # Check for and warn about zero coverage positions
    zero_cov_count = (df["total_coverage"] == 0).sum()
    if zero_cov_count > 0:
        logger.warning(
            f"Found {zero_cov_count:,} positions with zero total_coverage - allele frequencies will be NaN"
        )

    # Flag rows with positive coverage; used to exclude zero-coverage rows
    # from present-only averages/stds while keeping include-zeros behavior explicit.
    df["has_cov"] = df["total_coverage"] > 0

    # Compute squared coverage for standard deviation calculations (using moment method)
    df["total_coverage_sq"] = df["total_coverage"] ** 2

    # Compute allele frequencies for ATGC
    for allele in ALLELES:
        df[f"{allele}_freq"] = df[allele].astype(float) / df["total_coverage"].where(
            df["total_coverage"] > 0, np.nan
        )
        df[f"{allele}_freq_sq"] = df[f"{allele}_freq"] ** 2

    return df


def compute_overall_statistics(
    df: pd.DataFrame, n_samples: int, treat_absent_as_zero: bool
) -> pd.DataFrame:
    """
    Compute overall mean and standard deviation statistics across all samples.

    Args:
        df: Combined DataFrame with all sample data
        n_samples: Total number of samples
        treat_absent_as_zero: Whether to treat absent positions as zero coverage

    Returns:
        DataFrame with overall statistics per position.

    Notes:
        - Coverage means/stds use a denominator that depends on `treat_absent_as_zero`:
            * If True, denominator = total number of samples (include absent as zeros)
            * If False, denominator = number of samples with positive coverage (present-only)
        - Allele frequencies (mean_freq_*, std_freq_*) are always computed using only the
          present samples as the denominator to avoid downward bias when many samples are absent.
    """

    logger.info("Calculating overall mean/std per position")
    # Aggregate coverage and allele frequency moments
    overall = _aggregate_with_freqs(
        df, ["contig", "position"], present_col_name="n_present", join_how="outer"
    )

    # Ensure n_present has no NaNs before downstream use/casting
    # overall["n_present"] = overall["n_present"].fillna(0).astype(int)

    # Denominator depends on mode: all samples (include zeros) or present only
    if treat_absent_as_zero:
        denom_overall = float(n_samples)  # Use total samples if including zeros
    else:
        denom_overall = overall["n_present"].astype(float)  # Use present count

    # Compute overall mean and std using moments
    overall_mean, overall_std = _std_from_moments(
        overall["sum_cov"], overall["sumsq_cov"], denom_overall
    )

    # Build result DataFrame with overall statistics
    # n_present always reflects actual samples with positive coverage, regardless of mode
    tmp = overall.assign(
        total_coverage=overall["sum_cov"],  # Total coverage across all samples
        mean_coverage=overall_mean,
        std_coverage=overall_std,
        n_samples=n_samples,
        n_present=overall["n_present"].astype(int),
    )

    result = tmp.reset_index()[
        [
            "contig",
            "position",
            "total_coverage",
            "mean_coverage",
            "std_coverage",
            "n_present",
            "n_samples",
        ]
    ]

    # Compute per-allele (ATGC) overall frequency stats
    # IMPORTANT: allele frequencies should be averaged over present samples only,
    # regardless of include_zeros mode to avoid biasing toward zero.
    denom_freq_overall = overall["n_present"].astype(float)
    allele_block_df = _compute_allele_means_stds(
        overall, denom_freq_overall
    ).reset_index()
    result = result.merge(allele_block_df, on=["contig", "position"], how="outer")

    return result


def compute_grouped_statistics(
    df: pd.DataFrame,
    meta: pd.DataFrame,
    result: pd.DataFrame,
    treat_absent_as_zero: bool,
) -> pd.DataFrame:
    """
    Compute grouped statistics for samples with group and/or time information.

    Args:
        df: Combined DataFrame with all sample data
        meta: Metadata DataFrame
        result: DataFrame with overall statistics to extend
        treat_absent_as_zero: Whether to treat absent positions as zero coverage

    Returns:
        Extended DataFrame with grouped statistics.

    Notes:
        - Coverage means/stds per group use a denominator that depends on `treat_absent_as_zero`:
            * If True, denominator = size of the group from metadata
            * If False, denominator = number of group samples with positive coverage (present-only)
        - Allele frequencies per group/time are always computed using only present samples in that
          group/time combination, regardless of `treat_absent_as_zero`.
    """

    # Determine which grouping columns are present in metadata
    grouping_cols = [
        c for c in ("group", "time") if c in meta.columns and meta[c].nunique() > 1
    ]
    if not grouping_cols:
        logger.warning(
            "No grouping columns with multiple values found - skipping grouped statistics"
        )
        return result

    # Convert grouping columns to categorical for performance
    for c in grouping_cols:
        df[c] = df[c].astype("category")
        meta[c] = meta[c].astype("category")

    # Iterate over grouping levels: first 'group', then 'group'+'time' if both
    for i in range(len(grouping_cols)):
        current_groups = grouping_cols[: i + 1]
        logger.info(f"Grouped stats for: {current_groups}")

        # Group by contig, position, and group(s), aggregate sums and counts
        idx = ["contig", "position", *current_groups]
        grouped = _aggregate_with_freqs(
            df, idx, present_col_name="n_present_group", join_how="outer"
        )

        # Get total group sizes from metadata for denominators
        group_sizes = meta.groupby(current_groups, observed=True).size()
        group_sizes_aligned = group_sizes.reindex(grouped.index.droplevel([0, 1]))

        # Set denominators and n_present based on mode
        # n_present_used always reflects actual samples with positive coverage (consistent with overall stats)
        n_present_used = grouped["n_present_group"]

        if treat_absent_as_zero:
            # Use full group size for averages
            denom = group_sizes_aligned.astype(float)
        else:
            # Use present count for averages
            denom = grouped["n_present_group"].astype(float)

        # Calculate means and stds for each group using moments
        means, stds = _std_from_moments(grouped["sum_cov"], grouped["sumsq_cov"], denom)

        # Pivot means and stds to wide format for each (contig, position)
        wide_mean = means.unstack(current_groups)
        wide_std = stds.unstack(current_groups)
        wide_npr = n_present_used.unstack(current_groups)

        ns_series = pd.Series(
            group_sizes_aligned.values,
            index=grouped.index,
            name="n_samples_group",
        )
        wide_ns = ns_series.unstack(current_groups)

        # Optional NA handling block retained as comments for clarity
        # if treat_absent_as_zero:
        #     wide_mean = wide_mean.fillna(0.0)
        #     wide_std = wide_std.fillna(0.0)
        #     for col in wide_ns.columns:
        #         if col in group_sizes.index:
        #             group_size = group_sizes[col]
        #             wide_ns[col] = wide_ns[col].fillna(group_size)

        # n_present always shows actual positive-coverage samples, regardless of mode
        # wide_npr = wide_npr.fillna(0).astype(int)

        # Rename columns to descriptive names
        mean_cols = _flatten_cols("mean_coverage_", current_groups, wide_mean.columns)
        std_cols = _flatten_cols("std_coverage_", current_groups, wide_std.columns)
        npr_cols = _flatten_cols("n_present_", current_groups, wide_npr.columns)
        ns_cols = _flatten_cols("n_samples_group_", current_groups, wide_ns.columns)

        # Apply the new column names
        wide_mean.columns = mean_cols
        wide_std.columns = std_cols
        wide_npr.columns = npr_cols
        wide_ns.columns = ns_cols

        # Collect all blocks for efficient concatenation
        blocks = [wide_mean, wide_std, wide_npr, wide_ns]

        # Compute grouped per-allele frequency stats and add to blocks
        # IMPORTANT: allele frequencies should be averaged over present samples only
        denom_freq_grouped = grouped["n_present_group"].astype(float)
        allele_stats = _compute_allele_means_stds(grouped, denom_freq_grouped)
        for allele in ALLELES:
            mean_col = f"mean_freq_{allele}"
            std_col = f"std_freq_{allele}"
            wide_mean_freq_allele = allele_stats[mean_col].unstack(current_groups)
            wide_std_freq_allele = allele_stats[std_col].unstack(current_groups)
            wide_mean_freq_allele.columns = _flatten_cols(
                f"mean_freq_{allele}_", current_groups, wide_mean_freq_allele.columns
            )
            wide_std_freq_allele.columns = _flatten_cols(
                f"std_freq_{allele}_", current_groups, wide_std_freq_allele.columns
            )
            blocks.append(wide_mean_freq_allele)
            blocks.append(wide_std_freq_allele)

        # Combine all blocks efficiently with concat and merge back
        block = pd.concat(blocks, axis=1).reset_index()
        result = result.merge(block, on=["contig", "position"], how="left")

    return result


# --- helpers: small, focused utilities ---------------------------------------
def _std_from_moments(sum_, sumsq, denom) -> Tuple[pd.Series, pd.Series]:
    """Compute mean and standard deviation from first and second moments.

    This uses the population variance identity Var(X) = E[X^2] - (E[X])^2.
    Inputs can be pandas Series/Index-aligned arrays or numpy arrays; ``denom``
    may be a scalar or a vector aligned to ``sum_``/``sumsq``.

    To avoid tiny negative values due to floating-point error, the variance is
    clamped at zero via ``np.maximum(var, 0.0)`` before taking the square root.

    Parameters:
        sum_ (pd.Series | np.ndarray | float): Sum of values per key.
        sumsq (pd.Series | np.ndarray | float): Sum of squared values per key.
        denom (pd.Series | np.ndarray | float): Denominator (count) used for the
            mean/std calculation. Can be a scalar (e.g., total sample size) or a
            vector per key (e.g., n_present per position/group).

    Returns:
        Tuple[pd.Series, pd.Series]: A pair (mean, std) computed elementwise.

    Examples:
        >>> import pandas as pd
        >>> sum_ = pd.Series([10.0, 6.0])
        >>> sumsq = pd.Series([110.0, 20.0])
        >>> denom = pd.Series([2.0, 2.0])
        >>> mean, std = _std_from_moments(sum_, sumsq, denom)
        >>> round(mean.iloc[0], 3), round(std.iloc[0], 3)
        (5.0, 5.477)
        >>> round(mean.iloc[1], 3), round(std.iloc[1], 3)
        (3.0, 1.0)

        Edge case with zero denominator yields NaN:
        >>> mean2, std2 = _std_from_moments(pd.Series([0.0]), pd.Series([0.0]), 0.0)
        >>> float(mean2.iloc[0]), float(std2.iloc[0])
        (nan, nan)
    """
    mean = sum_ / denom
    var = (sumsq / denom) - mean**2
    return mean, np.sqrt(np.maximum(var, 0.0))


def _sanitize(x: object) -> str:
    """Sanitize labels for safe column naming.

    Replaces any character not in [A-Za-z0-9_.-] with an underscore. Input is
    first converted to ``str``. Useful for turning group/time labels into clean
    column suffixes.

    Parameters:
        x (object): Any value representing a label.

    Returns:
        str: A sanitized string with only [A-Za-z0-9_.-] and underscores.

    Examples:
        >>> _sanitize("Group A / time-1")
        'Group_A_time-1'
        >>> _sanitize(42)
        '42'
        >>> _sanitize("a:b(c)d")
        'a_b_c_d'
    """
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def _flatten_cols(prefix: str, current_groups: List[str], cols) -> List[str]:
    """Convert a MultiIndex-like list of group labels into flat column names.

    If an element of ``cols`` is a tuple (e.g., from unstacked MultiIndex), it is
    paired with the corresponding group names from ``current_groups`` to form
    parts like ``group_A`` and ``time_1``. If it is a scalar label, only the first
    entry in ``current_groups`` is used. All labels are sanitized by ``_sanitize``.

    Parameters:
        prefix (str): Column name prefix to start with, e.g., 'mean_coverage_'.
        current_groups (List[str]): Ordered list of group keys, e.g., ['group', 'time'].
        cols (Iterable[Any]): Labels from an unstacked index, either tuples or scalars.

    Returns:
        List[str]: Flattened column names of the form
            '<prefix>group_<g>_time_<t>' or '<prefix>group_<g>' when a single
            grouping dimension is present.

    Examples:
        Single grouping level:
        >>> _flatten_cols('mean_coverage_', ['group'], ['A', 'B'])
        ['mean_coverage_group_A', 'mean_coverage_group_B']

        Two grouping levels (tuples expected in cols):
        >>> _flatten_cols('std_coverage_', ['group', 'time'], [('A', '1'), ('B', '2')])
        ['std_coverage_group_A_time_1', 'std_coverage_group_B_time_2']
    """
    out = []
    for lab in cols:
        if isinstance(lab, tuple):
            parts = [f"{c}_{_sanitize(v)}" for c, v in zip(current_groups, lab)]
        else:
            parts = [f"{current_groups[0]}_{_sanitize(lab)}"]
        out.append(prefix + "_".join(parts))
    return out


def _aggregate_with_freqs(
    df: pd.DataFrame,
    idx: List[str],
    present_col_name: str,
    join_how: str = "outer",
) -> pd.DataFrame:
    """Aggregate coverage moments and allele-frequency moments on an index.

    Produces sum of coverage, sum of squared coverage, and count of present samples
    (based on ``has_cov``). Additionally joins sums of allele frequencies and their
    squares for each allele in ALLELES.

    Parameters:
        df: Long-form DataFrame containing 'total_coverage', 'total_coverage_sq',
            'has_cov', and per-allele frequency columns '<A>_freq' and
            '<A>_freq_sq' for A in ALLELES.
        idx: Grouping columns to aggregate by, e.g., ['contig', 'position'] or
            ['contig', 'position', 'group'].
        present_col_name: Name to use for the present-count column in the result
            (e.g., 'n_present' or 'n_present_group').
        join_how: How to join frequency aggregates onto coverage aggregates; kept
            as 'outer' by default to match existing behavior.

    Returns:
        Aggregated DataFrame indexed by ``idx`` with columns:
        ['sum_cov', 'sumsq_cov', present_col_name, '<A>_freq_sum', '<A>_freq_sumsq', ...]
    """
    cov_agg = (
        df.groupby(idx, observed=True)
        .agg(
            sum_cov=("total_coverage", "sum"),
            sumsq_cov=("total_coverage_sq", "sum"),
            present=("has_cov", "sum"),
        )
        .sort_index()
    )
    cov_agg = cov_agg.rename(columns={"present": present_col_name})

    # Build allele frequency aggregation dict
    agg_add = {}
    for allele in ALLELES:
        agg_add[f"{allele}_freq_sum"] = (f"{allele}_freq", "sum")
        agg_add[f"{allele}_freq_sumsq"] = (f"{allele}_freq_sq", "sum")

    freq_agg = df.groupby(idx, observed=True).agg(**agg_add).sort_index()
    return cov_agg.join(freq_agg, how=join_how)


def _compute_allele_means_stds(
    agg: pd.DataFrame, denom, alleles: List[str] = ALLELES
) -> pd.DataFrame:
    """Compute per-allele mean and std from aggregated frequency moments.

    Parameters:
        agg: Aggregated DataFrame containing '<A>_freq_sum' and '<A>_freq_sumsq'.
        denom: Denominator (scalar or aligned Series) used for moments.
        alleles: List of allele labels to compute stats for.

    Returns:
        DataFrame with columns 'mean_freq_<A>' and 'std_freq_<A>' for each allele,
        indexed like ``agg``.
    """
    out = pd.DataFrame(index=agg.index)
    for allele in alleles:
        mean_allele_freq, std_allele_freq = _std_from_moments(
            agg[f"{allele}_freq_sum"], agg[f"{allele}_freq_sumsq"], denom
        )
        out[f"mean_freq_{allele}"] = mean_allele_freq
        out[f"std_freq_{allele}"] = std_allele_freq
    return out


# --- main computation function ------------------------------------------------


def compute_mean_coverage(
    mag_metadata: Path, treat_absent_as_zero: bool = False
) -> pd.DataFrame:
    """
    Computes mean and standard deviation coverage statistics for a single MAG.

    This function orchestrates the computation of coverage statistics by:
    1. Loading and combining sample data from multiple profile files
    2. Computing allele frequencies for statistical analysis
    3. Calculating overall statistics across all samples
    4. Computing grouped statistics if group/time metadata is available

    Args:
        mag_metadata: Path to the MAG's metadata TSV file.
        treat_absent_as_zero: If True, positions not found in a sample are
            treated as having zero coverage. If False, the average is taken
            only over samples where the position is present.

    Returns:
        A pandas DataFrame containing the mean coverage and allele frequency statistics,
        with columns for contig, position, overall statistics, and group-specific statistics.
    """
    # Read and validate the MAG metadata file
    meta = read_mag_metadata(mag_metadata)
    n_samples = len(meta)
    logger.info(f"Aggregating coverage across {n_samples} samples from {mag_metadata}")

    # Load and combine all sample data
    df = load_and_combine_sample_data(meta)

    # Handle case where no data was found
    if df.empty:
        logger.warning("No coverage data found in any sample. Returning empty result.")
        return pd.DataFrame(
            columns=[
                "contig",
                "position",
                "total_coverage",
                "mean_coverage",
                "std_coverage",
                "n_present",
                "n_samples",
            ]
        )

    # Compute allele frequencies and prepare data for statistical calculations
    df = compute_allele_frequencies(df)
    # Log the mode for accounting of absent positions
    mode_desc = "include-zeros" if treat_absent_as_zero else "present-only"
    logger.info(f"Grouped mode: {mode_desc}")

    # Compute overall statistics across all samples
    result = compute_overall_statistics(df, n_samples, treat_absent_as_zero)

    # Compute grouped statistics if group/time information is available
    result = compute_grouped_statistics(df, meta, result, treat_absent_as_zero)

    # Finalize and return results
    base = [
        "contig",
        "position",
        "total_coverage",
        "mean_coverage",
        "std_coverage",
        "n_present",
        "n_samples",
    ]

    # Collect all group-related columns dynamically by prefix
    groupish = [
        col
        for col in result.columns
        if col.startswith(
            (
                "mean_coverage_",  # total coverage grouped means
                "std_coverage_",  # total coverage grouped stds
                "n_present_",
                "n_samples_group_",
                "mean_freq_",
                "std_freq_",
            )
        )
        and col not in {"mean_coverage", "std_coverage"}
    ]

    # Combine base and group columns, sort by contig and position, round to 6 decimals
    out = result[base + groupish].sort_values(["contig", "position"]).round(6)
    return out


def process_metadata_file_worker(
    metadata_file: str, output_dir: str, treat_absent_as_zero: bool
) -> str:
    """
    Worker function to process a single metadata file in parallel.

    This function is designed to be called by multiprocessing.Pool workers.
    It processes one MAG metadata file, computes coverage statistics, and saves the output.

    Args:
        metadata_file: Path to the metadata TSV file for a single MAG.
        output_dir: Directory where output files should be saved.
        treat_absent_as_zero: Whether to treat absent positions as zero coverage.

    """
    mag_id = os.path.basename(metadata_file).split("_metadata")[0]

    logger.info(f"Processing MAG: {mag_id}")

    # Compute mean coverage for this MAG
    df = compute_mean_coverage(
        mag_metadata=Path(metadata_file),
        treat_absent_as_zero=treat_absent_as_zero,
    )

    # Save output file
    out_file = Path(output_dir) / f"{mag_id}_mean_coverage.tsv"
    df.to_csv(out_file, sep="\t", index=False)

    logger.info(
        f"Saved mean coverage ({df.shape[0]:,} rows) for {mag_id} to {out_file}"
    )


def main():
    """
    Main entry point for the script.

    Parses command-line arguments, finds all *_metadata.tsv files in the input directory,
    and processes each MAG to compute mean coverage per position. Outputs are saved as
    TSV files in the specified output directory.

    This function uses multiprocessing to process multiple metadata files concurrently,
    utilizing all available CPU cores for improved performance when processing many files.
    """
    # Initialize logging for the module
    setup_logging()
    # Set up command-line argument parser
    args = argparse.ArgumentParser(
        description=(
            "Calculate mean per-position coverage across samples for MAGs "
            "using per-MAG metadata TSV files."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    args.add_argument(
        "--rootDir",
        required=True,
        type=str,
        help="Directory containing per-MAG *_metadata.tsv files.",
    )
    args.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="Directory to write output TSV files.",
    )
    args.add_argument(
        "--include_zeros",
        action="store_true",
        help=(
            "Include absent positions as zero coverage when calculating averages (default: average only over present samples)."
        ),
    )
    args.add_argument(
        "--cpus",
        type=int,
        default=multiprocessing.cpu_count(),
        help=(
            f"Number of CPU cores to use for parallel processing, all available cores)."
        ),
    )
    args = args.parse_args()

    treat_absent_as_zero = args.include_zeros

    os.makedirs(args.output_dir, exist_ok=True)
    metadata_files = glob(os.path.join(args.rootDir, "*_metadata.tsv"))
    if not metadata_files:
        logger.error(
            f"No *_metadata.tsv files found in input directory: {args.rootDir}"
        )
        return

    num_files = len(metadata_files)
    # Ensure we don't use more CPUs than files to process
    num_cpus = min(args.cpus, num_files)

    logger.info(f"Found {num_files} metadata file(s) to process")
    logger.info(f"Using {num_cpus} CPU core(s) for parallel processing")

    # Prepare task arguments for starmap
    tasks = [
        (metadata_file, args.output_dir, treat_absent_as_zero)
        for metadata_file in metadata_files
    ]

    # Process files using multiprocessing with starmap; exceptions will propagate
    with multiprocessing.Pool(processes=num_cpus) as pool:
        pool.starmap(process_metadata_file_worker, tasks)


if __name__ == "__main__":
    main()
