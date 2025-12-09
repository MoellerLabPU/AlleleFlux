#!/usr/bin/env python
"""
Generate preprocessing eligibility tables by aggregating status files from preprocessing steps.

This script reads JSON status files produced by preprocess_between_groups.py and
preprocess_within_group.py, then generates eligibility tables that indicate which
MAGs have sufficient positions remaining after filtering to proceed to downstream
statistical tests.
"""

import argparse
import json
import logging
import os
from glob import glob
from pathlib import Path

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def aggregate_between_groups_status(
    status_dir: str, min_positions: int = 1
) -> pd.DataFrame:
    """
    Aggregate status files from between-groups preprocessing.

    Parameters:
        status_dir: Directory containing *_preprocessing_status.json files
        min_positions: Minimum positions required for eligibility (overrides value in status file)

    Returns:
        DataFrame with columns:
            MAG_ID, input_rows, output_rows, total_positions,
            positions_eligible_unpaired, positions_eligible_paired,
            two_sample_unpaired_eligible, two_sample_paired_eligible
    """
    status_files = glob(os.path.join(status_dir, "*_preprocessing_status.json"))

    if not status_files:
        logger.warning(f"No status files found in {status_dir}")
        return pd.DataFrame(
            columns=[
                "MAG_ID",
                "input_rows",
                "output_rows",
                "total_positions",
                "positions_eligible_unpaired",
                "positions_eligible_paired",
                "two_sample_unpaired_eligible",
                "two_sample_paired_eligible",
            ]
        )

    records = []
    for status_file in status_files:
        with open(status_file, "r") as f:
            status = json.load(f)

        total_positions = status["total_positions"]
        positions_eligible_unpaired = status["positions_eligible_unpaired"]
        positions_eligible_paired = status["positions_eligible_paired"]

        # Use position-level eligibility with min_positions threshold
        # Unpaired tests (two_sample_unpaired, LMM) use unpaired eligibility
        # Paired tests (two_sample_paired, CMH) use paired eligibility
        eligible_unpaired = positions_eligible_unpaired >= min_positions
        eligible_paired = positions_eligible_paired >= min_positions

        records.append(
            {
                "MAG_ID": status["mag_id"],
                "input_rows": status["input_rows"],
                "output_rows": status["output_rows"],
                "total_positions": total_positions,
                "positions_eligible_unpaired": positions_eligible_unpaired,
                "positions_eligible_paired": positions_eligible_paired,
                # Unpaired tests (two_sample_unpaired, LMM) use unpaired eligibility
                "two_sample_unpaired_eligible": eligible_unpaired,
                # Paired tests (two_sample_paired, CMH) use paired eligibility
                "two_sample_paired_eligible": eligible_paired,
            }
        )

    df = pd.DataFrame(records)
    logger.info(
        f"Aggregated {len(records)} status files from between-groups preprocessing"
    )
    logger.info(
        f"Eligible MAGs (unpaired): {df['two_sample_unpaired_eligible'].sum()} / {len(df)}"
    )
    logger.info(
        f"Eligible MAGs (paired): {df['two_sample_paired_eligible'].sum()} / {len(df)}"
    )

    return df


def aggregate_within_groups_status(
    status_dir: str, groups: list, min_positions: int = 1
) -> pd.DataFrame:
    """
    Aggregate status files from within-groups preprocessing.

    Parameters:
        status_dir: Directory containing *_preprocessing_status.json files
        groups: List of group names (e.g., ["fat", "control"])
        min_positions: Minimum positions required for eligibility

    Returns:
        DataFrame with columns:
            MAG_ID, total_positions_{group1}, positions_eligible_{group1},
            single_sample_eligible_{group1}, ... (repeated for each group)
    """
    status_files = glob(os.path.join(status_dir, "*_preprocessing_status.json"))

    # Build column names
    base_columns = ["MAG_ID"]
    for group in groups:
        base_columns.extend(
            [
                f"total_positions_{group}",
                f"positions_eligible_{group}",
                f"single_sample_eligible_{group}",
            ]
        )

    if not status_files:
        logger.warning(f"No status files found in {status_dir}")
        return pd.DataFrame(columns=base_columns)

    # Group status files by MAG
    mag_statuses = {}
    for status_file in status_files:
        with open(status_file, "r") as f:
            status = json.load(f)

        mag_id = status["mag_id"]
        group = status["group"]
        total_positions = status["total_positions"]
        positions_eligible = status["positions_eligible"]

        # Use position-level eligibility with min_positions threshold
        eligible = positions_eligible >= min_positions

        if mag_id not in mag_statuses:
            mag_statuses[mag_id] = {"MAG_ID": mag_id}

        # Store total_positions, positions_eligible count and eligibility for each group
        mag_statuses[mag_id][f"total_positions_{group}"] = total_positions
        mag_statuses[mag_id][f"positions_eligible_{group}"] = positions_eligible
        # All within-group tests (single_sample, lmm_across_time, cmh_across_time) share eligibility
        mag_statuses[mag_id][f"single_sample_eligible_{group}"] = eligible

    df = pd.DataFrame(list(mag_statuses.values()))

    # Ensure all expected columns exist (fill with 0/False for missing groups)
    for col in base_columns:
        if col not in df.columns:
            if "positions" in col:
                df[col] = 0
            else:
                df[col] = False

    # Reorder columns
    df = df[base_columns]

    logger.info(
        f"Aggregated status files for {len(mag_statuses)} MAGs from within-groups preprocessing"
    )
    for group in groups:
        col = f"single_sample_eligible_{group}"
        if col in df.columns:
            logger.info(f"  {group}: {df[col].sum()} / {len(df)} MAGs eligible")

    return df


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Generate preprocessing eligibility tables from status files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--status_dir",
        required=True,
        help="Directory containing preprocessing status JSON files",
        type=str,
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Output path for eligibility table TSV",
        type=str,
    )
    parser.add_argument(
        "--preprocess_type",
        required=True,
        choices=["between_groups", "within_groups"],
        help="Type of preprocessing to aggregate",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        default=[],
        help="Group names (required for within_groups type)",
    )
    parser.add_argument(
        "--min_positions",
        type=int,
        default=1,
        help="Minimum number of positions required after filtering for eligibility",
    )

    args = parser.parse_args()

    if args.preprocess_type == "between_groups":
        df = aggregate_between_groups_status(args.status_dir, args.min_positions)
    else:
        if not args.groups:
            raise ValueError(
                "--groups is required for within_groups preprocessing type"
            )
        df = aggregate_within_groups_status(
            args.status_dir, args.groups, args.min_positions
        )

    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    df.to_csv(args.output_file, sep="\t", index=False)
    logger.info(f"Eligibility table written to {args.output_file}")


if __name__ == "__main__":
    main()
