#!/usr/bin/env python3
"""Copy MAG profile files for selected samples based on metadata filtering.

This script copies MAG profile files from a source profiles directory to a
destination, filtering samples based on group and optionally timepoint criteria
from the metadata file.
"""

import argparse
import logging
import os
import shutil
import sys
from glob import glob
from pathlib import Path
from typing import List, Optional

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def read_and_filter_metadata(
    metadata_file: str,
    groups: List[str],
    timepoints: Optional[List[str]] = None,
) -> List[str]:
    """Read metadata file and filter samples by group and optionally timepoint.

    Args:
        metadata_file: Path to the metadata TSV file.
        groups: List of group names to include.
        timepoints: Optional list of timepoints to include. If None, no time filtering.

    Returns:
        List of unique sample IDs matching the filter criteria.
    """
    metadata = pd.read_csv(metadata_file, sep="\t")

    # Determine required columns based on filtering needs
    required_columns = {"sample_id", "group"}
    if timepoints:
        required_columns.add("time")

    if not required_columns.issubset(metadata.columns):
        missing = required_columns - set(metadata.columns)
        logger.error(f"Missing required columns in metadata: {', '.join(missing)}")
        sys.exit(1)

    # Apply group filter (always required)
    filtered = metadata[metadata["group"].isin(groups)]

    # Apply timepoint filter if specified
    if timepoints:
        filtered = filtered[filtered["time"].isin(timepoints)]

    if filtered.empty:
        logger.error("No samples match the specified criteria.")
        sys.exit(1)

    sample_ids = filtered["sample_id"].unique().tolist()
    return sample_ids


def get_mag_ids_from_sample(
    profiles_dir: Path, sample_id: str
) -> List[str]:
    """Get all available MAG IDs from a sample's profile directory.

    Args:
        profiles_dir: Root profiles directory.
        sample_id: Sample ID (subdirectory name).

    Returns:
        List of MAG IDs found in the sample directory.
    """
    sample_dir = profiles_dir / sample_id
    if not sample_dir.is_dir():
        return []

    # Pattern: {sample_id}_{mag_id}_profiled.tsv.gz
    profile_files = list(sample_dir.glob(f"{sample_id}_*_profiled.tsv.gz"))
    mag_ids = []
    for f in profile_files:
        # Extract MAG ID from filename: {sample_id}_{mag_id}_profiled.tsv.gz
        name = f.name
        prefix = f"{sample_id}_"
        suffix = "_profiled.tsv.gz"
        if name.startswith(prefix) and name.endswith(suffix):
            mag_id = name[len(prefix) : -len(suffix)]
            mag_ids.append(mag_id)
    return mag_ids


def copy_mag_profiles(
    root_dir: Path,
    sample_ids: List[str],
    mag_ids: Optional[List[str]],
    destination: Path,
) -> None:
    """Copy MAG profile files from source to destination.

    Args:
        root_dir: Root directory containing 'profiles' subdirectory.
        sample_ids: List of sample IDs to process.
        mag_ids: List of MAG IDs to copy. If None, copy all MAGs.
        destination: Destination root directory (profiles subdir will be created).
    """
    src_profiles_dir = root_dir / "profiles"
    dest_profiles_dir = destination / "profiles"

    if not src_profiles_dir.is_dir():
        logger.error(f"Source profiles directory not found: {src_profiles_dir}")
        sys.exit(1)

    dest_profiles_dir.mkdir(parents=True, exist_ok=True)

    total_files_copied = 0
    total_files_not_found = 0
    samples_processed = 0

    for sample_id in sample_ids:
        sample_src_dir = src_profiles_dir / sample_id

        if not sample_src_dir.is_dir():
            logger.warning(
                f"Sample directory '{sample_src_dir}' does not exist. Skipping."
            )
            continue

        samples_processed += 1

        # Determine which MAGs to copy for this sample
        if mag_ids is None:
            # Copy all MAGs for this sample
            mags_to_copy = get_mag_ids_from_sample(src_profiles_dir, sample_id)
            if not mags_to_copy:
                logger.warning(f"No MAG profiles found in '{sample_src_dir}'.")
                continue
        else:
            mags_to_copy = mag_ids

        files_copied = 0
        files_not_found = 0

        for mag_id in mags_to_copy:
            file_name = f"{sample_id}_{mag_id}_profiled.tsv.gz"
            src_file = sample_src_dir / file_name

            if src_file.is_file():
                # Create destination sample directory
                dest_sample_dir = dest_profiles_dir / sample_id
                dest_sample_dir.mkdir(parents=True, exist_ok=True)

                shutil.copy2(src_file, dest_sample_dir / file_name)
                files_copied += 1
            else:
                logger.debug(f"File not found: {src_file}")
                files_not_found += 1

        total_files_copied += files_copied
        total_files_not_found += files_not_found

        if files_copied > 0 or files_not_found > 0:
            logger.debug(
                f"Sample '{sample_id}': copied {files_copied}, not found {files_not_found}"
            )

    logger.info("Copying Summary:")
    logger.info(f"  Samples processed: {samples_processed}/{len(sample_ids)}")
    logger.info(f"  Files copied: {total_files_copied}")
    if total_files_not_found > 0:
        logger.info(f"  Files not found: {total_files_not_found}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description=(
            "Copy MAG profile files for selected samples based on metadata filtering. "
            "Copies from {rootDir}/profiles/{sample_id}/ to {destination}/profiles/{sample_id}/."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rootDir",
        type=str,
        required=True,
        help="Root directory containing 'profiles' subdirectory with sample subdirectories.",
    )
    parser.add_argument(
        "--metadata_file",
        type=str,
        required=True,
        help="Path to the metadata TSV file (must have 'sample_id' and 'group' columns).",
    )
    parser.add_argument(
        "--destination",
        type=str,
        required=True,
        help="Destination root directory (profiles subdirectory will be created).",
    )
    parser.add_argument(
        "--groups",
        type=str,
        nargs="+",
        required=True,
        help="Groups to include (e.g., --groups control treatment).",
    )
    parser.add_argument(
        "--timepoints",
        type=str,
        nargs="*",
        default=None,
        help="Timepoints to filter by (optional). If not provided, all timepoints included.",
    )
    parser.add_argument(
        "--mag_ids",
        type=str,
        nargs="*",
        default=None,
        help="MAG IDs to copy. If not provided, all MAGs for each sample are copied.",
    )
    args = parser.parse_args()

    # Filter samples based on metadata
    sample_ids = read_and_filter_metadata(
        args.metadata_file,
        groups=args.groups,
        timepoints=args.timepoints if args.timepoints else None,
    )
    logger.info(f"Filtered {len(sample_ids)} samples based on metadata criteria.")
    logger.info(f"  Groups: {args.groups}")
    if args.timepoints:
        logger.info(f"  Timepoints: {args.timepoints}")

    if args.mag_ids:
        logger.info(f"Copying {len(args.mag_ids)} specified MAG(s).")
    else:
        logger.info("Copying all available MAGs for each sample.")

    copy_mag_profiles(
        root_dir=Path(args.rootDir),
        sample_ids=sample_ids,
        mag_ids=args.mag_ids,
        destination=Path(args.destination),
    )


if __name__ == "__main__":
    main()
