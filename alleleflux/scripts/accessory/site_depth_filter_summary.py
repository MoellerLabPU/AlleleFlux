#!/usr/bin/env python3
"""Run-level summary of the per-site depth filter.

Every ``alleleflux-cache-allele-freq`` job writes a sidecar next to its parquet,
``{mag}_{timepoint}_allele_frequency_site_depth_filter.tsv``, with one row per
sample: the mean coverage the sample was judged against and how many positions
the filter kept or removed.  This command concatenates every sidecar under the
run's ``allele_freq_cache/`` directory into one table and rolls it up per
MAG x timepoint, so "how much did the filter take out, and where" is one file.

Example: sidecars for MAG_A at t0 (samples s1, s2) and at t1 (s1) ->
``<output>`` with 3 rows (the sidecars stacked) and ``<output stem>_by_mag.tsv``
with 2 rows::

    mag_id  timepoint  n_samples  mean_of_mean_coverage  mean_of_mean_coverage_after_filter  sites_removed_high_total  sites_covered_total  fraction_removed_high
    MAG_A   t0         2          25.0                   22.0                                40                        3000                 0.0133
    MAG_A   t1         1          40.0                   40.0                                0                         1000                 0.0

``fraction_removed_high`` is the share of covered positions (summed over the
samples) that the ceiling removed: the number to look at when deciding whether
``site_ratio`` is biting too hard on a MAG.
"""

import argparse
import glob
import logging
import os

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

SIDECAR_GLOB = os.path.join("*", "*_site_depth_filter.tsv")


def main(argv=None):
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Concatenate the allele-frequency cache's site-depth-filter sidecars "
        "into one run-level table and a per-MAG x timepoint rollup.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--cache_dir",
        required=True,
        help="The run's allele_freq_cache/ directory (one subdirectory per timepoint).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Run-level TSV; the rollup is written next to it with the suffix _by_mag.tsv.",
    )
    args = parser.parse_args(argv)

    files = sorted(glob.glob(os.path.join(args.cache_dir, SIDECAR_GLOB)))
    if not files:
        # A run with the filter wired in always has sidecars; none means the
        # wrong directory or a cache built before the filter existed.
        raise FileNotFoundError(
            f"no *_site_depth_filter.tsv sidecars under {args.cache_dir}"
        )
    # sample_id / mag_id / timepoint are labels: read them as str so numeric
    # looking labels (group "20", timepoint "5") never become integers.
    run = pd.concat(
        [
            pd.read_csv(
                f, sep="\t", dtype={"sample_id": str, "mag_id": str, "timepoint": str}
            )
            for f in files
        ],
        ignore_index=True,
    )
    run.to_csv(args.output, sep="\t", index=False)

    # Per MAG x timepoint: the ceiling's bite as a fraction of covered positions.
    roll = run.groupby(["mag_id", "timepoint"], as_index=False).agg(
        n_samples=("sample_id", "nunique"),
        mean_of_mean_coverage=("mean_coverage", "mean"),
        mean_of_mean_coverage_after_filter=("mean_coverage_after_filter", "mean"),
        sites_removed_high_total=("n_removed_high", "sum"),
        sites_covered_total=("n_sites_covered", "sum"),
    )
    roll["fraction_removed_high"] = (
        roll["sites_removed_high_total"] / roll["sites_covered_total"]
    )
    roll_path = os.path.splitext(args.output)[0] + "_by_mag.tsv"
    roll.to_csv(roll_path, sep="\t", index=False)
    logger.info(
        f"{len(run):,} sample rows from {len(files)} sidecars -> {args.output}; "
        f"{len(roll)} MAG x timepoint rows -> {roll_path}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
