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
    """
    Perform a paired t-test on two related samples of scores.

    Parameters:
    group1 (array-like): The first set of scores.
    group2 (array-like): The second set of scores.

    Returns:
    float: The p-value of the t-test.
    """
    t_stat, t_p = stats.ttest_rel(
        group1,
        group2,
        nan_policy="raise",
        alternative="two-sided",
    )
    return t_p


def run_wilcoxon(group1, group2):
    """
    Perform the Wilcoxon signed-rank test on two related samples.

    Parameters:
    group1 (array-like): The first set of observations.
    group2 (array-like): The second set of observations.

    Returns:
    float: The p-value of the test.

    Raises:
    ValueError: If the input arrays contain NaN values.
    """
    w_stat, w_p = stats.wilcoxon(
        group1,
        group2,
        nan_policy="raise",
        alternative="two-sided",
    )
    return w_p


def paired_test_for_nucleotide(values, test_type):
    """
    Conducts a paired statistical test (t-test or Wilcoxon test) on a list of nucleotide values.

    Parameters:
    values (list or array-like): A list or array of nucleotide values to be tested.
    test_type (str): The type of test to perform. Options are "t-test", "wilcoxon", "both", or "either".

    Returns:
    tuple: A tuple containing the p-values of the tests. If test_type is "t-test", returns (t_p, None).
           If test_type is "wilcoxon", returns (None, w_p). If test_type is "both" or "either", returns (t_p, w_p).

    Notes:
    - If the number of values is less than or equal to 3, the function returns (1, None) for "t-test" and (None, 1) for "wilcoxon".
    - If the number of values is odd, the value closest to the median is removed before conducting the test.
    - If the values in both groups are identical, the function returns (1, None) for "t-test" and (None, 1) for "wilcoxon".
    - The function uses `run_tTest` and `run_wilcoxon` to perform the actual statistical tests.
    """
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


def process_group(group_df, test_type, data_type="longitudinal"):
    """
    Processes a group of nucleotide data and performs paired tests for each nucleotide.

    Parameters:
        group_df (pd.DataFrame): A DataFrame containing nucleotide data.
        test_type (str): The type of test to perform on the nucleotide values.
        data_type (str): The type of data ('longitudinal' or 'single').

    Returns:
        dict: A dictionary where keys are nucleotides and values are lists of p-values
              resulting from the paired tests.

    Raises:
        ValueError: If any NA values are found in the required columns or if any NaN p-values
                    are encountered in the results.
    """
    pvalues = {}
    for nuc in NUCLEOTIDES:
        # Column name varies depending on data_type
        if data_type == "longitudinal":
            col = f"{nuc}_diff_mean"
        else:  # single
            col = f"{nuc}"

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
    """
    Processes a genomic site to determine if it should be removed based on statistical tests.
    Only remove the site only if ALL 4 nuleotides have nuc_remove as True.
    If any one nucleotide has nuc_remove as False, then remove_site is False.

    Parameters:
        args (tuple): A tuple containing:
            - (contig, position) (tuple): The contig and position of the site.
            - group_df (DataFrame): A DataFrame containing the group data for the site.
            - alpha (float): The significance level for the statistical tests.
            - test_type (str): The type of statistical test to perform. Can be "t-test", "wilcoxon", "both", or "either".
            - data_type (str): The type of data ('longitudinal' or 'single').

    Returns:
        tuple: A tuple containing:
            - (contig, position) (tuple): The contig and position of the site.
            - remove_site (bool): A boolean indicating whether the site should be removed (True) or not (False).
    """
    (contig, position), group_df, alpha, test_type, data_type = args
    pvals = process_group(group_df, test_type=test_type, data_type=data_type)
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

        if not nuc_remove:
            remove_site = False
            break
    return ((contig, position), remove_site)


def filter_sites_parallel(grouped, alpha, test_type, cpus, data_type="longitudinal"):
    """
    Filters sites in parallel using multiprocessing.

    Parameters:
        grouped (iterable): An iterable of grouped data, where each element is a tuple
                            ((contig, position), group_df).
        alpha (float): The significance level for the statistical test.
        test_type (str): The type of statistical test to perform.
        cpus (int): The number of CPU cores to use for parallel processing.
        data_type (str): The type of data ('longitudinal' or 'single').

    Returns:
        list: A list of sites to remove, where each site is represented by its (contig, position).
    """
    groups = list(grouped)  # list of ((contig, position), group_df)
    args_list = [
        ((contig, position), group_df, alpha, test_type, data_type)
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
        "--data_type",
        help="Is the data from a single timepoint or from a time series (longitudinal)?",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
    )
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
    df = pd.read_csv(args.mean_changes_fPath, sep="\t", dtype={"gene_id": str})

    # Group data by "contig" and "position" (including groups with NA keys).
    grouped = df.groupby(["contig", "position"], dropna=False)

    logging.info(
        f"Determining sites to remove based on paired tests using {args.cpus} cpus. Test type is set to {args.test_type}. Data type is {args.data_type}"
    )
    sites_to_remove = filter_sites_parallel(
        grouped, args.alpha, args.test_type, args.cpus, args.data_type
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
