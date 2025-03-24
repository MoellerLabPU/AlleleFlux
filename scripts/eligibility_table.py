#!/usr/bin/env python3
import argparse
import glob
import logging
import os
import sys

import numpy as np
import pandas as pd


def generate_test_eligibility_table(qc_dir, min_sample_num, data_type):
    # Get all QC files in the directory that match *_QC.tsv
    file_list = glob.glob(os.path.join(qc_dir, "*_QC.tsv"))
    if not file_list:
        raise ValueError("No QC files found in the specified directory.")

    # Read all files into a list of DataFrames
    df_list = []
    for f in file_list:
        logging.info(f"Reading QC file: {f}")
        df = pd.read_csv(f, sep="\t")
        df_list.append(df)

    # Concatenate all QC DataFrames
    qc_df = pd.concat(df_list, ignore_index=True)

    # Group by MAG_ID and group to get one set of counts per combination.
    # We assume the values in replicates_per_group and paired_replicates_per_group are constant for each group.
    summary = qc_df.groupby(["MAG_ID", "group"], as_index=False).agg(
        {"replicates_per_group": "first", "paired_replicates_per_group": "first"}
    )

    # Pivot the table so that each row is a unique combination of MAG_ID and the groups become separate columns.
    pivot = summary.pivot(index=["MAG_ID"], columns="group")
    # This creates a MultiIndex for columns; flatten it:
    pivot.columns = ["_".join([str(col[0]), str(col[1])]) for col in pivot.columns]
    pivot = pivot.reset_index()

    # Get the list of groups from the original summary (should be exactly 2)
    groups = summary["group"].unique()
    if len(groups) != 2:
        raise ValueError(f"Expected exactly 2 groups, found: {groups}")
    g1, g2 = groups

    # Check unpaired test eligibility: True if replicates_per_group for both groups > min_sample_num.
    pivot["unpaired_test_eligible"] = (
        pivot[f"replicates_per_group_{g1}"] >= min_sample_num
    ) & (pivot[f"replicates_per_group_{g2}"] >= min_sample_num)

    # Check paired test eligibility: True if paired_replicates_per_group for both groups >= min_sample_num.
    pivot["paired_test_eligible"] = (
        pivot[f"paired_replicates_per_group_{g1}"] >= min_sample_num
    ) & (pivot[f"paired_replicates_per_group_{g2}"] >= min_sample_num)

    # For single-sample test eligibility, compute per group based on data_type:
    if data_type == "longitudinal":
        pivot[f"single_sample_eligible_{g1}"] = (
            pivot[f"replicates_per_group_{g1}"] >= min_sample_num
        )
        pivot[f"single_sample_eligible_{g2}"] = (
            pivot[f"replicates_per_group_{g2}"] >= min_sample_num
        )
    else:
        # Set single sample eligibility columns to NaN for single timepointdata
        pivot[f"single_sample_eligible_{g1}"] = np.nan
        pivot[f"single_sample_eligible_{g2}"] = np.nan

    return pivot


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Generate test eligibility summary table from QC files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--qc_dir", required=True, help="Directory containing *_QC.tsv files", type=str
    )
    parser.add_argument(
        "--min_sample_num",
        type=int,
        default=4,
        help="Minimum number of replicates per group required",
    )
    parser.add_argument(
        "--output_file",
        required=True,
        help="Output file for eligibility summary table",
        type=str,
    )
    
    parser.add_argument(
        "--data_type",
        help="Is the data from a single timepoint or from a time series (longitudinal)",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )

    args = parser.parse_args()

    eligibility_table = generate_test_eligibility_table(
        args.qc_dir, args.min_sample_num, args.data_type
    )
    eligibility_table.to_csv(args.output_file, sep="\t", index=False)
    logging.info(f"Eligibility summary table written to {args.output_file}")


if __name__ == "__main__":
    main()
