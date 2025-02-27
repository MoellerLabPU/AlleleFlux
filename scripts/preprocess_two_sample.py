#!/usr/bin/env python
import argparse
import logging
import time
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import supress_warning
from scipy import stats
from tqdm import tqdm

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]


def run_tTest(group1, group2):
    t_stat, t_p = stats.ttest_rel(
        group1,
        group2,
        nan_policy="raise",
        alternative="two-sided",
    )
    return t_p


def run_wilcoxon(group1, group2):
    w_stat, w_p = stats.wilcoxon(
        group1,
        group2,
        nan_policy="raise",
        alternative="two-sided",
    )
    return w_p


def paired_test_for_nucleotide(values, test_type):
    # Assign p-values of 1 if there are less than 3 values. List of 3 will be converted to 2.
    # paired test using 2 value for t-test generally gives NaN and 1 for wilcoxon
    if len(values) <= 3:
        if test_type == "t-test":
            return (1, None)
        elif test_type == "wilcoxon":
            return (None, 1)
        else:
            return (1, 1)

    # Handle odd numbers: remove the value closest to the median.
    if len(values) % 2 != 0:
        median = np.median(values)
        distances = np.abs(values - median)
        # Find all indices with min distance
        min_indices = np.where(distances == distances.min())[0]
        # Remove first occurrence of closest to median
        values = np.delete(values, min_indices[0])

    # np.sort sorts in ascending order, so we use [::-1] to sort in descending order
    sorted_vals = np.sort(values)[::-1]

    # Compute the midpoint index of the sorted array.
    half = len(sorted_vals) // 2
    # First half of sorted values
    group1 = sorted_vals[:half]
    # Second half of sorted values
    group2 = sorted_vals[half:]

    # Check for identical values
    d = group1 - group2
    if np.all(d == 0):
        if test_type == "t-test":
            # Paired t-test outputs NaN if both groups have identical values
            return (1, None)
        elif test_type == "wilcoxon":
            # Wilcoxon gives an error if both groups have identical values
            return (None, 1)
        else:
            return (1, 1)

    # Conduct paired t-test.
    if test_type == "t-test":
        t_p = run_tTest(group1, group2)
        return (t_p, None)
    elif test_type == "wilcoxon":
        w_p = run_wilcoxon(group1, group2)
        return (None, w_p)
    elif test_type in ["both", "either"]:
        t_p = run_tTest(group1, group2)
        w_p = run_wilcoxon(group1, group2)
        return (t_p, w_p)


def process_group(group_df, test_type):
    pvalues = {}
    for nuc in NUCLEOTIDES:
        col = f"{nuc}_diff_mean"
        values = group_df[col].values
        if pd.isnull(values).any():
            raise ValueError(f"NA values found in column '{col}'")
        # pvalues[nuc] = paired_test_for_nucleotide(values, test_type=test_type)
        nuc_pvals = paired_test_for_nucleotide(values, test_type=test_type)
        # Check that none of the returned p-values are NaN (ignoring None values)
        for p in nuc_pvals:
            if p is not None and np.isnan(p):
                raise ValueError(
                    f"NaN p-value encountered for nucleotide {nuc} at site: "
                    f"{group_df[['contig', 'position']].iloc[0].to_dict()}"
                )
        pvalues[nuc] = nuc_pvals
    return pvalues


def process_site(args):
    (contig, position), group_df, alpha, test_type = args
    pvals = process_group(group_df, test_type=test_type)
    remove_site = True  # assume removal unless one nucleotide prevents it
    for nuc, (t_p, w_p) in pvals.items():
        if test_type == "t-test":
            nuc_remove = t_p >= alpha
        elif test_type == "wilcoxon":
            nuc_remove = w_p >= alpha
        elif test_type == "both":
            # Remove the site if both tests are not significant
            nuc_remove = (t_p >= alpha) and (w_p >= alpha)
        elif test_type == "either":
            # Remove if either test is non-significant.
            nuc_remove = (t_p >= alpha) or (w_p >= alpha)
        else:
            raise ValueError(f"Unknown test type: {test_type}")
        if not nuc_remove:
            remove_site = False
            break
    return ((contig, position), remove_site)


def filter_sites_parallel(grouped, alpha, test_type, cpus):
    groups = list(grouped)  # list of ((contig, position), group_df)
    args_list = [
        ((contig, position), group_df, alpha, test_type)
        for (contig, position), group_df in groups
    ]
    with Pool(processes=cpus) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(process_site, args_list),
                total=len(args_list),
                desc="Processing sites",
            )
        )
    sites_to_remove = [site for site, remove_flag in results if remove_flag]
    return sites_to_remove


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    parser = argparse.ArgumentParser(
        description="Filter sites using paired tests on sorted values from mean_changes_df.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--mean_changes_fPath",
        required=True,
        help="Path to mean changes dataframe",
        type=str,
    )

    parser.add_argument("--output_fPath", required=True, help="Path to output file")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance level")
    parser.add_argument(
        "--test_type",
        choices=["t-test", "wilcoxon", "both", "either"],
        default="t-test",
        help="""
        - "t-test":   remove site if t_p >= alpha \n
        - "wilcoxon": remove site if w_p >= alpha \n
        - "both":     remove if (t_p >= alpha) AND (w_p >= alpha)
                    (i.e., remove if both tests are non-significant)
        - "either":   remove if (t_p >= alpha) OR (w_p >= alpha)
                    (i.e., remove if at least one test is non-significant)
        """,
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=cpu_count(),
        help="Number of processors to use",
    )
    args = parser.parse_args()
    start_time = time.time()

    logging.info("Reading input file")
    df = pd.read_csv(args.mean_changes_fPath, sep="\t")

    # Group data by "contig" and "position" (including groups with NA keys).
    grouped = df.groupby(["contig", "position"], dropna=False)

    logging.info(
        f"Determining sites to remove based on paired tests using {args.cpus} cpus."
    )
    sites_to_remove = filter_sites_parallel(
        grouped, args.alpha, args.test_type, args.cpus
    )
    logging.info(f"Removing {len(sites_to_remove):,} sites")

    # Remove these sites from the DataFrame.
    removal_index = pd.MultiIndex.from_tuples(
        sites_to_remove, names=["contig", "position"]
    )
    df_filtered = (
        df.set_index(["contig", "position"])
        .drop(removal_index, errors="raise")
        .reset_index()
    )

    df_filtered.to_csv(args.output_fPath, sep="\t", index=False)
    logging.info(
        f"Retained {len(df_filtered):,} rows. Filtered data written to {args.output_fPath}"
    )
    end_time = time.time()
    logging.info(f"Total time taken: {end_time-start_time:.2f} seconds")


if __name__ == "__main__":
    main()
