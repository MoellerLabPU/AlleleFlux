#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import sys

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def read_and_filter_metadata(metadata_file):
    metadata = pd.read_csv(metadata_file, sep="\t")

    required_columns = {"sample", "time", "group"}
    if not required_columns.issubset(metadata.columns):
        missing = required_columns - set(metadata.columns)
        logger.error(f"Missing required columns in metadata: {', '.join(missing)}")
        sys.exit(1)

    # Apply filters
    filtered = metadata[
        (metadata["time"].isin(["pre", "end", "post"]))
        & (metadata["group"].isin(["control", "fat"]))
    ]

    if filtered.empty:
        logger.error("No samples match the specified criteria.")
        sys.exit(1)

    # Extract sample IDs
    sample_ids = filtered["sample"].unique()
    return sample_ids


def copy_mag_profiles_structured(rootDir, sample_ids, mag_ids, destination):
    # Ensure the destination root directory exists
    os.makedirs(destination, exist_ok=True)

    total_files_copied = 0
    total_files_not_found = 0

    for mag_id in mag_ids:
        logger.info(f"Processing MAG_ID: {mag_id}")
        files_copied = 0
        files_not_found = 0

        for sample_id in sample_ids:
            sorted_dir = f"{sample_id}.sorted"
            sorted_dir_path = os.path.join(rootDir, sorted_dir)

            if not os.path.isdir(sorted_dir_path):
                logger.warning(
                    f"Directory '{sorted_dir_path}' does not exist. Skipping sample '{sample_id}'."
                )
                files_not_found += 1
                continue

            # Construct the expected file name
            file_name = f"{sorted_dir}_{mag_id}_profiled.tsv.gz"
            file_path = os.path.join(sorted_dir_path, file_name)

            if os.path.isfile(file_path):
                # Create corresponding directory in destination
                dest_sorted_dir = os.path.join(destination, sorted_dir)
                os.makedirs(dest_sorted_dir, exist_ok=True)

                shutil.copy2(file_path, dest_sorted_dir)
                files_copied += 1
            else:
                logger.warning(f"File '{file_name}' not found in '{sorted_dir_path}'.")
                files_not_found += 1

        logger.info(f"Summary for MAG_ID '{mag_id}':")
        logger.info(f"  Total samples processed: {len(sample_ids)}")
        logger.info(f"  Files successfully copied: {files_copied}")
        logger.info(f"  Files not found or failed to copy: {files_not_found}")

        total_files_copied += files_copied
        total_files_not_found += files_not_found

    logger.info("Overall Copying Summary:")
    logger.info(f"Total MAG_IDs processed: {len(mag_ids)}")
    logger.info(f"Total samples processed: {len(sample_ids) * len(mag_ids)}")
    logger.info(f"Total files successfully copied: {total_files_copied}")
    logger.info(f"Total files not found or failed to copy: {total_files_not_found}")


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description=(
            "Copy MAG profile files for samples where "
            "time is 'end' or 'post' and group is 'control' or 'fat', "
            "replicating the directory structure in the destination."
        )
    )
    parser.add_argument(
        "--rootDir",
        type=str,
        required=True,
        help="Root directory containing .sorted subdirectories.",
    )
    parser.add_argument(
        "--metadata_file",
        type=str,
        required=True,
        help="Path to the metadata file.",
    )
    parser.add_argument(
        "--destination",
        type=str,
        required=True,
        help="Destination root directory where structured MAG files will be copied.",
    )
    parser.add_argument(
        "--mag_ids",
        type=str,
        nargs="+",
        required=True,
        help="MAG ID to filter files.",
    )
    args = parser.parse_args()

    sample_ids = read_and_filter_metadata(args.metadata_file)
    logger.info(f"Filtered {len(sample_ids)} samples based on metadata criteria.")

    copy_mag_profiles_structured(
        rootDir=args.rootDir,
        sample_ids=sample_ids,
        mag_ids=args.mag_ids,
        destination=args.destination,
    )


if __name__ == "__main__":
    main()
