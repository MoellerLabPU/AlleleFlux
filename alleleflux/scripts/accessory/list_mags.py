#!/usr/bin/env python3
"""
Generate a list of MAG IDs from a directory of metadata files.

This utility script extracts MAG identifiers from *_metadata.tsv files
and writes them to a text file (one per line). Useful for:
- Creating input for Snakemake workflows
- Batch processing scripts
- Generating MAG subsets for testing

Usage:
    alleleflux-list-mags --metadata_dir /path/to/metadata --output mag_list.txt
    alleleflux-list-mags --metadata_dir /path/to/metadata --output mag_list.txt --qc_dir /path/to/qc --min_samples 5
"""

import argparse
import logging
import os
from glob import glob
from pathlib import Path

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)


def get_mag_list(metadata_dir: str, qc_dir: str = None, min_samples: int = 1) -> list:
    """
    Extract MAG IDs from metadata directory with optional filtering.

    Args:
        metadata_dir: Directory containing *_metadata.tsv files
        qc_dir: Optional directory with QC files to filter by passing samples
        min_samples: Minimum number of samples (after QC) required to include MAG

    Returns:
        List of MAG IDs meeting criteria
    """
    metadata_files = glob(os.path.join(metadata_dir, "*_metadata.tsv"))

    if not metadata_files:
        logger.warning(f"No *_metadata.tsv files found in {metadata_dir}")
        return []

    logger.info(f"Found {len(metadata_files)} metadata file(s)")

    mag_ids = []
    for meta_file in metadata_files:
        mag_id = os.path.basename(meta_file).replace("_metadata.tsv", "")

        # Read metadata to count samples
        meta_df = pd.read_csv(meta_file, sep="\t")
        n_samples = len(meta_df)

        # If QC filtering is enabled, count only passing samples
        if qc_dir:
            qc_file = Path(qc_dir) / f"{mag_id}_QC.tsv"
            if qc_file.exists():
                qc_df = pd.read_csv(qc_file, sep="\t")
                n_passing = (qc_df["breadth_threshold_passed"] == True).sum()
                logger.debug(f"{mag_id}: {n_passing}/{n_samples} samples passed QC")
                n_samples = n_passing
            else:
                logger.warning(f"QC file not found for {mag_id}, using all samples")

        # Check if MAG meets minimum sample threshold
        if n_samples >= min_samples:
            mag_ids.append(mag_id)
            logger.debug(f"Including {mag_id} ({n_samples} samples)")
        else:
            logger.debug(f"Excluding {mag_id} ({n_samples} < {min_samples} samples)")

    logger.info(f"Selected {len(mag_ids)} MAG(s) meeting criteria")
    return sorted(mag_ids)


def main():
    """Main entry point for MAG list generation."""
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Generate list of MAG IDs from metadata directory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--metadata_dir",
        required=True,
        type=str,
        help="Directory containing *_metadata.tsv files",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        help="Output file path for MAG list (one ID per line)",
    )
    parser.add_argument(
        "--qc_dir",
        type=str,
        help="Directory with QC files to filter by passing samples only",
    )
    parser.add_argument(
        "--min_samples",
        type=int,
        default=1,
        help="Minimum number of samples (after QC) required to include MAG",
    )

    args = parser.parse_args()

    # Validate input directory
    if not os.path.isdir(args.metadata_dir):
        logger.error(f"Metadata directory does not exist: {args.metadata_dir}")
        return 1

    # Validate QC directory if provided
    if args.qc_dir and not os.path.isdir(args.qc_dir):
        logger.error(f"QC directory does not exist: {args.qc_dir}")
        return 1

    # Get MAG list
    mag_ids = get_mag_list(
        metadata_dir=args.metadata_dir,
        qc_dir=args.qc_dir,
        min_samples=args.min_samples,
    )

    if not mag_ids:
        logger.error("No MAGs found meeting criteria")
        return 1

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for mag_id in mag_ids:
            f.write(f"{mag_id}\n")

    logger.info(f"Wrote {len(mag_ids)} MAG ID(s) to {output_path}")
    return 0


if __name__ == "__main__":
    exit(main())
