import argparse
import logging
import os

import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging

# Constants
NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]

logger = logging.getLogger(__name__)


def preprocess_data(df, max_zero_count, output_dir, mag_id, group):
    """
    Preprocess data by removing positions where more than max_zero_count replicates
    have _diff_mean value zero for all nucleotides.

    Parameters:
    -----------
    df : pandas.DataFrame
        Input dataframe containing nucleotide frequency data
    max_zero_count : int
        Maximum allowed replicates per position where all nucleotide_diff_mean values are zero
    output_dir : str
        Directory to save the processed data
    mag_id : str
        MAG identifier
    group : str
        Group identifier

    Returns:
    --------
    pandas.DataFrame
        Processed dataframe with positions removed based on zero count threshold
    """
    logger.info(
        f"Preprocessing: Removing positions where more than {max_zero_count} replicates have _diff_mean value zero for all nucleotides."
    )

    # List the columns we want to check.
    diff_mean_cols = [f"{nuc}_diff_mean" for nuc in NUCLEOTIDES]

    # For each replicate, check if all nucleotide diff_mean columns are zero.
    df["all_zeros"] = (df[diff_mean_cols] == 0).all(axis=1)

    # Group by position (using contig and position) and count replicates with all zeros.
    # No need to include gene_id here as we are not using that again for the final dataframe as with the significance tests
    zeros_count = df.groupby(["contig", "position"], dropna=False)["all_zeros"].sum()

    # Identify positions where the number of replicates with all zeros exceeds the threshold.
    positions_to_remove = zeros_count[zeros_count > max_zero_count].index
    logger.info(
        f"Removing {len(positions_to_remove):,} positions based on zero count threshold."
    )

    # Remove these positions from the dataframe.
    df = df.set_index(["contig", "position"])
    df = df.drop(positions_to_remove).reset_index()

    # Cleanup temporary column
    df = df.drop(columns=["all_zeros"])

    df.to_csv(
        os.path.join(
            output_dir,
            f"{mag_id}_{group}_allele_frequency_changes_mean_zeros_processed.tsv.gz",
        ),
        index=False,
        sep="\t",
        compression="gzip",
    )


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Script to preprocess within group data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--df_fPath",
        required=True,
        help="Path to mean changes dataframe",
        type=str,
    )
    parser.add_argument(
        "--max_zero_count",
        type=int,
        default=4,
        help="Maximum allowed replicates per position where all nucleotide_diff_mean values are zero. Positions where the count of such replicates exceeds this threshold are removed.",
    )
    parser.add_argument(
        "--mag_id",
        help="MAG ID to process",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--group",
        help="Group to perform tests on",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    input_df = pd.read_csv(args.df_fPath, sep="\t", dtype={"gene_id": str})

    # Filter the dataframe to only include positions from the specified group.
    input_df = input_df[input_df["group"] == args.group]
    logger.info(
        f"Data filtered for group '{args.group}', resulting in {input_df.shape[0]:,} rows."
    )

    os.makedirs(args.output_dir, exist_ok=True)

    preprocess_data(
        df=input_df,
        max_zero_count=args.max_zero_count,
        output_dir=args.output_dir,
        mag_id=args.mag_id,
        group=args.group,
    )


if __name__ == "__main__":
    main()
