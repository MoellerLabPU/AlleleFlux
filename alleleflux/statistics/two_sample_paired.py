import argparse
import gc
import logging
import os
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

import alleleflux.utilities.supress_warning as supress_warning

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]


def perform_paired_tests(func_paired, grouped_df, cpus, num_tests, output_dir, mag_id):
    start_time = time.time()

    with Pool(processes=cpus) as pool:
        results_iter = pool.imap_unordered(func_paired, grouped_df)
        records = []
        for result in tqdm(
            results_iter, desc="Performing paired significance tests", total=num_tests
        ):
            name_tuple, p_values, num_pairs, notes = result
            contig, gene_id, position = name_tuple
            record = {
                "contig": contig,
                "gene_id": gene_id,
                "position": position,
                **p_values,
                "num_pairs": num_pairs,
                "notes": notes,
            }
            records.append(record)

    end_time = time.time()
    logging.info(f"Paired tests performed in {end_time - start_time:.2f} seconds")
    test_results = pd.DataFrame(records)
    # Identify p-value columns
    p_value_columns = [col for col in test_results.columns if "_p_value" in col]

    # Remove rows where all p-value columns are NaN
    test_results.dropna(subset=p_value_columns, how="all", inplace=True)

    logging.info(
        f"Saving 2 sample paired significance results for MAG {mag_id} to {output_dir}"
    )
    test_results.to_csv(
        os.path.join(output_dir, f"{mag_id}_two_sample_paired.tsv.gz"),
        index=False,
        sep="\t",
        compression="gzip",
    )


def run_paired_tests(args, group_1, group_2, min_sample_num, data_type="longitudinal"):
    name_tuple, grouped_df = args
    # Separate the data into two groups
    group1 = grouped_df[grouped_df["group"] == group_1]
    group2 = grouped_df[grouped_df["group"] == group_2]

    # Merge the two groups on 'replicate_id', 'contig', and 'position'
    merged_data = pd.merge(
        group1,
        group2,
        on=["replicate", "contig", "position"],
        suffixes=("_group1", "_group2"),
        how="inner",
    )

    # Initialize p_values with NaN
    p_values = {}
    notes = ""
    for nucleotide in NUCLEOTIDES:
        p_values[f"{nucleotide}_p_value_paired_tTest"] = np.nan
        p_values[f"{nucleotide}_p_value_Wilcoxon"] = np.nan

    # Number of pairs
    num_pairs = merged_data.shape[0]

    # Only perform the tests if there are at least min_sample_num pairs
    if num_pairs >= min_sample_num:
        for nucleotide in NUCLEOTIDES:
            # Use either diff_mean columns (longitudinal) or direct values (single)
            if data_type == "longitudinal":
                data1 = merged_data[f"{nucleotide}_diff_mean_group1"]
                data2 = merged_data[f"{nucleotide}_diff_mean_group2"]
            else:  # Single data
                data1 = merged_data[f"{nucleotide}_group1"]
                data2 = merged_data[f"{nucleotide}_group2"]

            # Check for identical values
            d = data1 - data2
            if np.all(d == 0):
                p_values[f"{nucleotide}_p_value_paired_tTest"] = (
                    1.0  # Paired t-test outputs NaN if both groups have identical values
                )
                p_values[f"{nucleotide}_p_value_Wilcoxon"] = (
                    1.0  # Wilcoxon gives an error if both groups have identical values
                )
                notes += f"{nucleotide}: identical values in both groups, p-value for both tests is set to 1; "
            else:
                # Perform paired t-test
                res_ttest = stats.ttest_rel(
                    data1,
                    data2,
                    nan_policy="raise",
                    alternative="two-sided",
                )
                p_values[f"{nucleotide}_p_value_paired_tTest"] = res_ttest.pvalue

                res_wilcoxon = stats.wilcoxon(
                    data1,
                    data2,
                    alternative="two-sided",
                    nan_policy="raise",
                )
                p_values[f"{nucleotide}_p_value_Wilcoxon"] = res_wilcoxon.pvalue

    return (name_tuple, p_values, num_pairs, notes)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Run two sample paired test for a MAG across different samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_df",
        required=True,
        help="Path to mean changes dataframe or allele frequency dataframe",
        type=str,
    )
    parser.add_argument(
        "--min_sample_num",
        type=int,
        default=4,
        help="Minimum number of replicates per group required",
    )

    parser.add_argument(
        "--mag_id",
        help="MAG ID to process",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--cpus",
        help=f"Number of processors to use.",
        default=cpu_count(),
        type=int,
    )

    parser.add_argument(
        "--output_dir",
        help="Path to output directory",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--data_type",
        help="Type of data to analyze: longitudinal or single",
        type=str,
        choices=["longitudinal", "single"],
        default="longitudinal",
    )

    args = parser.parse_args()

    input_file = args.input_df
    logging.info(f"Loading data from {input_file}")

    # Load the input data
    logging.info("Reading input dataframe..")
    input_df = pd.read_csv(input_file, sep="\t", dtype={"gene_id": str})

    # Get unique groups
    groups = input_df["group"].unique()
    if len(groups) != 2:
        raise ValueError(
            f"Expected exactly 2 groups for 2-sample tests tests, but found {len(groups)} groups: {groups}. Exiting...."
        )

    group_1, group_2 = groups

    # Group the data
    logging.info("Grouping data by contig, gene_id and position")
    grouped_df = input_df.groupby(["contig", "gene_id", "position"], dropna=False)

    num_tests = len(grouped_df)
    os.makedirs(args.output_dir, exist_ok=True)

    logging.info(
        f"Performing {num_tests:,} paired tests (paired by replicate) between {group_1} and {group_2} using {args.cpus} cores."
    )
    func_paired = partial(
        run_paired_tests,
        group_1=group_1,
        group_2=group_2,
        min_sample_num=args.min_sample_num,
        data_type=args.data_type,
    )
    perform_paired_tests(
        func_paired, grouped_df, args.cpus, num_tests, args.output_dir, args.mag_id
    )


if __name__ == "__main__":
    main()
