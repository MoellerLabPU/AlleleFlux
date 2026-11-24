#!/usr/bin/env python3
"""Per-(MAG, timepoint) allele-frequency cache writer.  [Stage 1 of 2]

Role in the two-stage pipeline
-------------------------------
This script is **Stage 1**.  It runs ONCE per (MAG, timepoint) regardless of
how many (timepoint_combination, group_combination) pairs include that
timepoint.

    alleleflux-cache-allele-freq
        --magID       MRGM_1079
        --qc_file     {OUTDIR}/QC/QC_5mo/MRGM_1079_QC.tsv
        --timepoint   5mo
        --data_type   longitudinal
        --output_path {OUTDIR}/allele_freq_cache/5mo/MRGM_1079_5mo_allele_frequency.parquet
        --cpus        16

Output consumed by
------------------
- ``allele_freq.py`` (Stage 2): reads exactly two cache files per combination
  to compute per-subject allele-frequency diffs.
- ``CMH.py``, ``LMM.py`` (across-time / between-groups tests): read one or
  two cache files directly when the longitudinal per-sample data is needed
  (replaces the old ``_allele_frequency_longitudinal.tsv.gz``).

Why this exists
---------------
Before this refactor, ``allele_freq.py`` was invoked once per
(timepoint_combination, group_combination).  For a config with timepoints
``[pre, post, end]`` and two combinations ``pre_post`` and ``pre_end``, the
``pre`` profile files were read, frequencies computed, and ``pre`` rows written
to disk **twice** — once inside each combination's longitudinal TSV.

This script eliminates that duplication: profile reads and frequency
computations for ``pre`` happen once; the resulting Parquet file is read by
both downstream combinations.

Shared helpers (from ``_allele_freq_common.py``)
-------------------------------------------------
- ``load_qc_results``        — load QC TSV and filter to passing samples
- ``init_worker``            — broadcast metadata dict to worker processes
- ``process_mag_files``      — read one profile TSV + compute frequencies (worker fn)
- ``build_metadata_from_qc`` — convert QC rows into pool-friendly structures

Note on QC input
-----------------
The Snakemake rule ``allele_freq_cache`` passes exactly ONE QC file per cache
job: the genome's QC table of the requested timepoint
(``QC/QC_{timepoint}/{mag}_QC.tsv``, one row per sample of that timepoint).

A genome with no usable sample (``--allow_empty``)
--------------------------------------------------
The workflow fixes the list of cache files before any QC table is read: under
the cohort baseline scope it asks for every tested genome at every witness
timepoint.  A genome can have no QC-passing sample at one of them (its QC
table there is header-only, or every sample failed).  With ``--allow_empty``
that job writes a cache with the normal columns and zero rows plus a
header-only sidecar, logs a warning and exits 0; the readers see a timepoint
with no sample.  Without the flag it is an error, as it is for a call by hand.
A sample that passed QC but whose profile cannot be read is an error either
way.
"""

import argparse
import logging
import re
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
    build_metadata_from_qc,
    empty_cache_table,
    init_worker,
    load_qc_results,
    process_mag_files,
)
from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.site_depth_filter import (
    add_site_depth_filter_args,
    site_depth_filter_from_args,
)

logger = logging.getLogger(__name__)

# Columns of the sidecar written next to every cache: one row per sample, what
# the site-depth filter kept and removed.  Spelled once so the sidecar of a
# zero-row cache (header only) has exactly the columns of a populated one.
SIDECAR_COLUMNS = [
    "sample_id",
    "mag_id",
    "timepoint",
    "mean_coverage",
    "mean_coverage_after_filter",
    "n_sites_covered",
    "n_removed_low",
    "n_removed_high",
    "n_kept",
]


def sidecar_path(output_path):
    """The sidecar's path for a cache path: same folder, ``.parquet`` replaced.

    Example: ``.../MAG_7_pre_allele_frequency.parquet`` ->
    ``.../MAG_7_pre_allele_frequency_site_depth_filter.tsv``.
    """
    return output_path.with_name(
        re.sub(r"\.parquet$", "", output_path.name) + "_site_depth_filter.tsv"
    )


def write_empty_cache(output_path, data_type):
    """Write a zero-row cache and its header-only sidecar (``--allow_empty``).

    The cache has the columns and types of a populated one
    (``empty_cache_table``), so every reader opens it like any other file and
    finds no sample in it.  Written to a temporary name and renamed, like a
    populated cache, so a killed job never leaves a file that looks complete.

    Parameters
    ----------
    output_path : pathlib.Path
        The cache file to write; its folder must exist.
    data_type : str
        ``"longitudinal"`` or ``"single"`` (decides whether ``time`` is a column).

    Returns
    -------
    pathlib.Path
        The sidecar that was written next to the cache.
    """
    tmp_path = output_path.with_suffix(".parquet.tmp")
    pq.write_table(
        empty_cache_table(data_type),
        tmp_path,
        compression="snappy",
        use_dictionary=True,
    )
    tmp_path.rename(output_path)
    sidecar = sidecar_path(output_path)
    pd.DataFrame(columns=SIDECAR_COLUMNS).to_csv(sidecar, sep="\t", index=False)
    return sidecar


def filter_qc_to_timepoint(qc_path, mag_id, timepoint):
    """Load one genome's QC table and keep the passing samples of one timepoint.

    The Snakemake rule supplies exactly one QC path per cache job: the
    genome's table of this timepoint, so every row already belongs to it.  The
    filter on ``time`` stays as a guard for a call by hand with a table that
    holds several timepoints.

    Parameters
    ----------
    qc_path : str
        Path to the per-MAG QC TSV (``QC/QC_{timepoint}/{mag}_QC.tsv``).
    mag_id : str
        Used for logging and error messages.
    timepoint : str
        Timepoint label to retain (matched against the ``time`` column after
        casting both sides to ``str`` to avoid int/str mismatches).

    Returns
    -------
    pd.DataFrame
        QC rows for samples at this timepoint that passed coverage/breadth QC.
        MAY BE EMPTY (a header-only table, or no sample passed): the caller
        decides whether that is an error (``--allow_empty``).

    Raises
    ------
    ValueError
        If the QC file has no ``coverage_threshold_passed`` or no ``time``
        column.
    """
    qc_df = load_qc_results(qc_path, mag_id)

    if "time" not in qc_df.columns:
        raise ValueError(
            f"QC file for MAG {mag_id} has no 'time' column; "
            "cannot filter to a single timepoint. This rule expects longitudinal QC data."
        )

    # Cast to str on both sides to handle numeric timepoints (e.g. 5 vs "5").
    qc_df["time"] = qc_df["time"].astype(str)
    qc_df = qc_df[qc_df["time"] == str(timepoint)]

    logger.info(
        f"Found {len(qc_df)} samples for MAG {mag_id} at timepoint '{timepoint}'."
    )
    return qc_df


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description=(
            "Stage 1: Compute per-sample allele frequencies for a MAG at a single "
            "timepoint and write to a Parquet cache file."
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
        "--qc_file",
        required=True,
        type=str,
        help=(
            "Path to the genome's QC TSV of this timepoint "
            "(QC/QC_{timepoint}/{mag}_QC.tsv).  One table: a cache is one "
            "genome at one timepoint."
        ),
    )
    parser.add_argument(
        "--timepoint",
        required=False,
        default=None,
        type=str,
        help=(
            "Timepoint label to extract from the QC file.  Required for "
            "'longitudinal' data_type (filters the QC file's 'time' column).  "
            "Not used for 'single' data_type — the QC file covers one timepoint "
            "and no row-level filtering is needed."
        ),
    )
    parser.add_argument(
        "--data_type",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
        help=(
            "Data type.  For 'longitudinal', rows in the QC file are filtered "
            "to --timepoint.  For 'single', all QC-passing samples are loaded."
        ),
    )
    parser.add_argument(
        "--output_path",
        required=True,
        type=str,
        help=(
            "Destination Parquet file.  Consumed by allele_freq.py (Stage 2) "
            "and by across-time statistical tests (CMH, LMM)."
        ),
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=cpu_count(),
        help="Worker processes for parallel profile loading.",
    )
    parser.add_argument(
        "--allow_empty",
        action="store_true",
        help=(
            "When the genome has NO QC-passing sample at this timepoint (a "
            "header-only QC table included), write a cache with the normal "
            "columns and zero rows plus a header-only sidecar, warn and exit 0. "
            "Without the flag that case is an error.  The workflow passes it "
            "because its list of cache files is fixed before QC is read.  A "
            "sample that passed QC but whose profile cannot be read stays an "
            "error."
        ),
    )
    add_site_depth_filter_args(parser)

    args = parser.parse_args()
    start_time = time.time()
    mag_id = args.magID

    # ------------------------------------------------------------------
    # 1. Determine which samples to load
    # ------------------------------------------------------------------
    if args.data_type == "longitudinal":
        if args.timepoint is None:
            parser.error("--timepoint is required when --data_type is 'longitudinal'.")
        # The Snakemake rule passes the genome's QC table of this timepoint.
        # Keep its passing rows at the requested timepoint.
        qc_df = filter_qc_to_timepoint(args.qc_file, mag_id, args.timepoint)
    else:
        # Single data: the QC file covers exactly one timepoint — load all passing
        # samples without any timepoint-based row filtering.
        qc_df = load_qc_results(args.qc_file, mag_id)

    output_path = Path(args.output_path)
    tp_label = args.timepoint if args.data_type == "longitudinal" else "N/A (single)"

    # No QC-passing sample at all.  An error for a call by hand; under
    # --allow_empty (the workflow) a zero-row cache, so the job still writes
    # the two files the workflow was promised.  This is decided HERE, from the
    # QC table alone, before any profile is opened: a profile that cannot be
    # read further down is a real failure and never ends up as an empty cache.
    if qc_df.empty:
        message = (
            f"zero QC-passing samples for MAG {mag_id} at timepoint {tp_label} "
            f"in {args.qc_file}"
        )
        if not args.allow_empty:
            raise ValueError(
                f"{message}; nothing to cache (pass --allow_empty to write a "
                f"zero-row cache instead)"
            )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        sidecar = write_empty_cache(output_path, args.data_type)
        logger.warning(
            f"{message}: empty cache written to {output_path} "
            f"(header-only sidecar {sidecar})"
        )
        return

    # ------------------------------------------------------------------
    # 2. Build worker inputs
    # ------------------------------------------------------------------
    # build_metadata_from_qc returns:
    #   metadata_dict: {sample_id -> metadata row} — broadcast to workers
    #   base_tuples:   [(sample_id, file_path), ...] — one entry per sample
    metadata_dict, base_tuples = build_metadata_from_qc(qc_df)

    # Append mag_id to each tuple to match the (sample_id, filepath, mag_id)
    # signature expected by process_mag_files.
    sample_file_tuples = [(sid, fpath, mag_id) for sid, fpath in base_tuples]

    # ------------------------------------------------------------------
    # 3. Load profiles in parallel and stream to Parquet
    # ------------------------------------------------------------------
    # Each worker:
    #   - reads one profile TSV (the expensive step) with int32 / string dtypes
    #   - attaches metadata (group, subjectID, replicate, time) as categoricals
    #   - computes per-base nucleotide frequencies as float32
    #
    # Rather than collect every per-sample frame in the parent and then
    # concat + write (which doubles parent-process memory at the moment of
    # concat), we stream each frame straight into a single Parquet file via
    # pyarrow.parquet.ParquetWriter.  Each per-sample table becomes one row
    # group; pyarrow handles per-row-group dictionaries for the categorical
    # columns automatically.
    n_proc = min(args.cpus, len(sample_file_tuples))
    logger.info(
        f"Computing allele frequencies for {len(sample_file_tuples)} samples "
        f"of MAG {mag_id} (timepoint={tp_label}) using {n_proc} processes."
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Write to a temp file first, then rename atomically.  If the job is
    # killed mid-write (SLURM SIGTERM/SIGKILL) or hits an OOM, the partial
    # file stays at the .tmp path and Snakemake correctly re-runs the job
    # because the final output is still missing.
    tmp_path = output_path.with_suffix(".parquet.tmp")

    writer = None
    total_rows = 0
    filter_summaries = []
    site_filter = site_depth_filter_from_args(args)
    logger.info(f"site-depth filter: {site_filter}")
    try:
        with Pool(
            processes=n_proc,
            initializer=init_worker,  # broadcasts metadata_dict + DATA_TYPE + the site filter
            initargs=(metadata_dict, args.data_type, site_filter),
        ) as pool:
            for df in pool.imap_unordered(process_mag_files, sample_file_tuples):
                # None means the worker logged an error for that sample;
                # skip it and continue with the rest.
                if df is None:
                    continue
                # One summary row per sample for the sidecar written below.
                filter_summaries.append(df.attrs["site_depth_filter"])

                # Convert to Arrow.  Categorical pandas columns map to Arrow
                # dictionary type; pyarrow then writes them as dictionary
                # pages in the row group.  preserve_index=False avoids an
                # extra __index_level_0__ column.
                table = pa.Table.from_pandas(df, preserve_index=False)

                if writer is None:
                    # Lock the schema on the first table.  All subsequent
                    # tables must match this schema (dtype + column order);
                    # process_mag_files produces them deterministically so
                    # this holds in practice.
                    writer = pq.ParquetWriter(
                        tmp_path,
                        table.schema,
                        compression="snappy",
                        use_dictionary=True,
                    )

                writer.write_table(table)
                total_rows += table.num_rows
                # Drop the per-sample frame ASAP — the next iteration will
                # build the next one and we don't want N frames live at once.
                del df, table
    finally:
        if writer is not None:
            writer.close()

    if total_rows == 0:
        logger.error(
            f"No sample profiles could be loaded for MAG {mag_id} "
            f"at timepoint {args.timepoint}."
        )
        # Remove the empty/missing-schema .tmp so the next attempt starts
        # from a clean slate.
        if tmp_path.exists():
            tmp_path.unlink()
        sys.exit(1)

    logger.info(f"Wrote {total_rows:,} rows to {output_path}")
    tmp_path.rename(output_path)

    # Sidecar next to the cache: what the site-depth filter removed, per
    # sample.  The run-level table is built from every sidecar by
    # alleleflux-site-filter-summary.
    sidecar = sidecar_path(output_path)
    summary_df = pd.DataFrame(filter_summaries)
    summary_df.insert(
        0, "timepoint", tp_label if args.data_type == "longitudinal" else ""
    )
    summary_df = summary_df[SIDECAR_COLUMNS]
    summary_df.to_csv(sidecar, sep="\t", index=False)
    logger.info(
        f"site-depth filter: {summary_df['n_removed_high'].sum():,} high-depth and "
        f"{summary_df['n_removed_low'].sum():,} low-depth sites removed over "
        f"{len(summary_df)} samples -> {sidecar}"
    )

    logger.info(f"Total time: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
