#!/usr/bin/env python3
"""Per-(MAG, timepoint) allele-frequency cache writer.  [Stage 1 of 2]

Role in the two-stage pipeline
-------------------------------
This script is **Stage 1**.  It runs ONCE per (MAG, timepoint) regardless of
how many (timepoint_combination, group_combination) pairs include that
timepoint.

    alleleflux-cache-allele-freq
        --magID       MRGM_1079
        --qc_files    QC_5mo_10mo-1D_AL/MRGM_1079_QC.tsv
        --timepoint   5mo
        --data_type   longitudinal
        --output_path allele_freq_cache/MRGM_1079_1D_AL_5mo_allele_frequency.parquet
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
The Snakemake rule ``compute_allele_freq_per_timepoint`` passes exactly ONE
canonical QC file per cache job (the QC from the first timepoint_combination
in config order that contains the requested timepoint).  ``filter_qc_to_timepoint``
handles this single-file case; there is no need to union across combinations.
"""

import argparse
import logging
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd

from alleleflux.scripts.analysis._allele_freq_common import (
    build_metadata_from_qc,
    init_worker,
    load_qc_results,
    process_mag_files,
)
from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def filter_qc_to_timepoint(qc_path, mag_id, timepoint):
    """Load the canonical QC file and filter to a single timepoint.

    The Snakemake rule always supplies exactly one QC path per cache job —
    the canonical QC file for this (gr_combo, timepoint) pair (the first
    timepoint_combination in config order that contains this timepoint).
    This function loads it, filters rows to ``timepoint``, and returns the
    result ready for ``build_metadata_from_qc``.

    Parameters
    ----------
    qc_path : str
        Path to the canonical per-MAG QC TSV for this (gr_combo, timepoint).
    mag_id : str
        Used for logging and error messages.
    timepoint : str
        Timepoint label to retain (matched against the ``time`` column after
        casting both sides to ``str`` to avoid int/str mismatches).

    Returns
    -------
    pd.DataFrame
        QC rows for samples at this timepoint that passed coverage/breadth QC.

    Raises
    ------
    ValueError
        If the QC file has no ``time`` column, or if no samples remain after
        filtering to the requested timepoint.
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

    if qc_df.empty:
        raise ValueError(
            f"No QC-passing samples for MAG {mag_id} at timepoint '{timepoint}'."
        )

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
        "--qc_files",
        required=True,
        nargs="+",
        type=str,
        help=(
            "Path(s) to per-MAG QC TSV file(s).  Exactly one file is expected "
            "for longitudinal data (the canonical QC file for this gr_combo and "
            "timepoint); one file for single data.  The nargs='+' form is kept "
            "so the Snakemake rule can expand the list naturally."
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

    args = parser.parse_args()
    start_time = time.time()
    mag_id = args.magID

    # ------------------------------------------------------------------
    # 1. Determine which samples to load
    # ------------------------------------------------------------------
    if args.data_type == "longitudinal":
        if args.timepoint is None:
            parser.error("--timepoint is required when --data_type is 'longitudinal'.")
        # The Snakemake rule passes one canonical QC file per (gr_combo, timepoint).
        # Filter its rows to the single requested timepoint.
        qc_df = filter_qc_to_timepoint(args.qc_files[0], mag_id, args.timepoint)
    else:
        # Single data: the QC file covers exactly one timepoint — load all passing
        # samples without any timepoint-based row filtering.
        qc_df = load_qc_results(args.qc_files[0], mag_id)

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
    # 3. Load profiles in parallel
    # ------------------------------------------------------------------
    # Each worker:
    #   - reads one profile TSV (the expensive step)
    #   - attaches metadata (group, subjectID, replicate, time)
    #   - computes per-base nucleotide frequencies (A/T/G/C ÷ total_coverage)
    # Returns a list of per-sample DataFrames.
    n_proc = min(args.cpus, len(sample_file_tuples))
    tp_label = args.timepoint if args.data_type == "longitudinal" else "N/A (single)"
    logger.info(
        f"Computing allele frequencies for {len(sample_file_tuples)} samples "
        f"of MAG {mag_id} (timepoint={tp_label}) using {n_proc} processes."
    )

    with Pool(
        processes=n_proc,
        initializer=init_worker,         # broadcasts metadata_dict + DATA_TYPE
        initargs=(metadata_dict, args.data_type),
    ) as pool:
        data_list = list(pool.imap_unordered(process_mag_files, sample_file_tuples))

    # Filter out None returns (failed samples — already logged as errors).
    data_list = [df for df in data_list if df is not None]
    if not data_list:
        logger.error(
            f"No sample profiles could be loaded for MAG {mag_id} "
            f"at timepoint {args.timepoint}."
        )
        sys.exit(1)

    # ------------------------------------------------------------------
    # 4. Write Parquet cache
    # ------------------------------------------------------------------
    # Concatenate all per-sample frames into one DataFrame.  The resulting
    # file holds every sample at this timepoint for this MAG — it is the
    # "cache" that Stage 2 (allele_freq.py) and statistical-test scripts read
    # instead of the old per-combination longitudinal TSV.
    combined = pd.concat(data_list, ignore_index=True)

    output_path = Path(args.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Writing {len(combined):,} rows to {output_path}")
    # Parquet with Snappy compression:
    #   - column-oriented → fast for downstream queries that select a subset of columns
    #   - Snappy is fast to decompress (better wall-clock than gzip for large reads)
    #   - typically 30–50% smaller on disk than gzipped TSV for this data
    combined.to_parquet(
        output_path,
        engine="pyarrow",
        compression="snappy",
        index=False,
    )

    logger.info(f"Total time: {time.time() - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
