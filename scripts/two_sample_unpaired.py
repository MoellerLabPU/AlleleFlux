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

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]


def perform_unpaired_tests(
    func_unpaired, grouped_df, group_1, group_2, cpus, num_tests, output_dir, mag_id
):
    start_time = time.time()

    with Pool(processes=cpus) as pool:
        results_iter = pool.imap_unordered(func_unpaired, grouped_df)
        records = []
        for result in tqdm(
            results_iter, desc="Performing significance tests", total=num_tests
        ):
            name_tuple, p_values, num_samples_group1, num_samples_group2, notes = result
            contig, gene_id, position = name_tuple
            record = {
                "contig": contig,
                "gene_id": gene_id,
                "position": position,
                **p_values,
                f"num_samples_{group_1}": num_samples_group1,
                f"num_samples_{group_2}": num_samples_group2,
                "notes": notes,
            }
            records.append(record)

    end_time = time.time()
    logging.info(
        f"2 sample unpaired tests performed in {end_time - start_time:.2f} seconds"
    )
    test_results = pd.DataFrame(records)
    # Identify p-value columns
    p_value_columns = [col for col in test_results.columns if "_p_value" in col]

    # Remove rows where all p-value columns are NaN
    test_results.dropna(subset=p_value_columns, how="all", inplace=True)

    logging.info(
        f"Saving 2 sample unpaired significance results for MAG {mag_id} to {output_dir}"
    )
    test_results.to_csv(
        os.path.join(output_dir, f"{mag_id}_two_sample_unpaired.tsv.gz"),
        index=False,
        sep="\t",
        compression="gzip",
    )


def run_unpaired_tests(
    args, group_1, group_2, min_sample_num, data_type="longitudinal"
):
    name_tuple, grouped_df = args
    # Separate the data into two groups
    group1 = grouped_df[grouped_df["group"] == group_1]
    group2 = grouped_df[grouped_df["group"] == group_2]

    # Initialize p_values with NaN
    p_values = {}
    notes = ""
    for nucleotide in NUCLEOTIDES:
        p_values[f"{nucleotide}_p_value_tTest"] = np.nan
        p_values[f"{nucleotide}_p_value_MannWhitney"] = np.nan

    # Counts of samples in each group
    num_samples_group1 = group1.shape[0]
    num_samples_group2 = group2.shape[0]

    # Only perform the t-test if both groups have at least min_sample_num data points
    if num_samples_group1 >= min_sample_num and num_samples_group2 >= min_sample_num:
        for nucleotide in NUCLEOTIDES:
            # Select appropriate column based on data type
            if data_type == "longitudinal":
                nuc_col = f"{nucleotide}_diff_mean"
            else:  # single data type
                nuc_col = nucleotide

            mean1 = np.mean(group1[nuc_col])
            mean2 = np.mean(group2[nuc_col])
            # ddof of 1 is used as we are calculating sample variance. "0" is used for population variance
            var1 = np.var(group1[nuc_col], ddof=1)
            var2 = np.var(group2[nuc_col], ddof=1)
            # In t-test if both groups have identical values, the p-value is NaN, eg. 5,5,5,5 and 5,5,5,5,5,5,5,5
            # I'm setting it 1 for consistency with Mann Whitney
            if var1 == 0 and var2 == 0 and mean1 == mean2:
                p_values[f"{nucleotide}_p_value_tTest"] = 1
                notes += f"{nucleotide}: identical values in both groups, p-value for t-test is set to 1; "
            else:
                # Perform t-test
                res_para = stats.ttest_ind(
                    group1[nuc_col],
                    group2[nuc_col],
                    equal_var=False,
                    nan_policy="raise",
                    alternative="two-sided",
                )
                p_values[f"{nucleotide}_p_value_tTest"] = res_para.pvalue

            # Perform Mann-Whitney U test
            res_non_para = stats.mannwhitneyu(
                group1[nuc_col],
                group2[nuc_col],
                alternative="two-sided",
                nan_policy="raise",  # No NaN values should be present. But if present a ValueError will be raised
            )
            p_values[f"{nucleotide}_p_value_MannWhitney"] = res_non_para.pvalue

    return (name_tuple, p_values, num_samples_group1, num_samples_group2, notes)


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Run two sample unpaired test for a MAG across different samples.",
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

    # Load the input data depending on the data type
    if args.data_type == "longitudinal":
        logging.info("Processing longitudinal data")
        input_df = pd.read_csv(input_file, sep="\t", dtype={"gene_id": str})
    else:
        logging.info("Processing single data")
        input_df = pd.read_csv(input_file, sep="\t", dtype={"gene_id": str})

    # Get unique groups
    groups = input_df["group"].unique()
    if len(groups) != 2:
        raise ValueError(
            f"Expected exactly 2 groups for 2-sample tests, but found {len(groups)} groups: {groups}. Exiting...."
        )

    group_1, group_2 = groups

    # Group the data
    logging.info("Grouping data by contig, gene_id and position")
    grouped_df = input_df.groupby(["contig", "gene_id", "position"], dropna=False)

    num_tests = len(grouped_df)

    logging.info(
        f"Performing {num_tests:,} unpaired tests between {group_1} and {group_2} using {args.cpus} cores."
    )
    os.makedirs(args.output_dir, exist_ok=True)

    func_unpaired = partial(
        run_unpaired_tests,
        group_1=group_1,
        group_2=group_2,
        min_sample_num=args.min_sample_num,
        data_type=args.data_type,
    )
    perform_unpaired_tests(
        func_unpaired,
        grouped_df,
        group_1,
        group_2,
        args.cpus,
        num_tests,
        args.output_dir,
        args.mag_id,
    )


if __name__ == "__main__":
    main()
