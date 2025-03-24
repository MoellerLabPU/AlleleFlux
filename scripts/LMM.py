#!/usr/bin/env python
import argparse
import logging
import warnings
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# import supress_warning
from tqdm import tqdm


def has_perfect_separation(sub_df, freq_col):
    """
    Determine if there is perfect separation in the dataset.

    Perfect separation occurs when a variable (or set of variables) can perfectly predict
    the outcome, which can cause issues in statistical models.

    This function checks if:
    1. All groups have zero variance (all values within a group are identical)
    2. There are at least two groups with different mean values

    Parameters
    ----------
    sub_df : pandas.DataFrame
        DataFrame containing the data to check for perfect separation.
        Must have a 'group' column for grouping.

    freq_col : str
        Name of the column to analyze for perfect separation within groups.

    Returns
    -------
    bool
        True if perfect separation is detected, False otherwise.
    """

    # Group by "group" and compute mean and variance
    group_stats = sub_df.groupby("group")[freq_col].agg(["mean", "var"])

    # Check:
    # 1. All groups have zero variance (var == 0)
    # 2. At least two groups have different means
    return (group_stats["var"] == 0).all() and (group_stats["mean"].nunique() > 1)


def detect_column_suffix(sub_df):
    # Identify if the columns are mean changes or raw frequencies
    base_freq_cols = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]
    freq_cols = []
    # Identify if the columns are mean changes or raw frequencies
    for col in base_freq_cols:
        if f"{col}_diff_mean" in sub_df.columns:
            freq_cols.append(f"{col}_diff_mean")
        elif col in sub_df.columns:
            freq_cols.append(col)
        else:
            raise ValueError(
                f"Neither '{col}' nor '{col}_diff_mean' found in dataframe columns."
            )
    return freq_cols


def run_model(args):
    contig, gene_id, position, sub_df = args

    position_result = {
        "contig": contig,
        "gene_id": gene_id,
        "position": position,
        "n_samples": len(sub_df),
        "notes": "",
    }

    def handle_case(nucleotide, note, p_value=np.nan, coef=np.nan, t_value=np.nan):
        position_result[f"{nucleotide}_coef"] = coef
        position_result[f"{nucleotide}_p_value"] = p_value
        position_result[f"{nucleotide}_t_value (z-score)"] = t_value
        position_result["notes"] += note

    freq_cols = detect_column_suffix(sub_df)

    for freq_col in freq_cols:
        nucleotide = freq_col.split("_")[0]
        warnings_list = []
        freq_values = sub_df[freq_col]

        with warnings.catch_warnings(record=True) as w:
            # https://docs.python.org/3/library/warnings.html#the-warnings-filter
            warnings.simplefilter("default")

            try:
                # Attempt model fitting
                model = smf.mixedlm(
                    f"{freq_col} ~ group", sub_df, groups=sub_df["replicate"]
                )
                result = model.fit()

            except np.linalg.LinAlgError as e:
                # Capture warnings first before continuing
                warnings_list = [str(warn.message) for warn in w]
                # These checks are done after model fitting to capture the warnings
                if has_perfect_separation(sub_df, freq_col):
                    handle_case(
                        nucleotide,
                        # P-value is tends towards 0, but LMM gives singular matrix error due to perfect separation
                        note=f"{nucleotide}: Singular matrix error, perfect separation, p-value set to 0; ",
                        p_value=0,
                    )
                    logging.info(
                        f"{contig}, {position}, {nucleotide}: Singular matrix error, perfect separation, p-value set to 0; ",
                    )
                elif freq_values.nunique() == 1:
                    # Set p-value as 1 if there are no unique values for the nucleotide
                    handle_case(
                        nucleotide,
                        note=f"{nucleotide}: Singular matrix error, constant values, p-value set to 1; ",
                        p_value=1,
                    )
                    logging.info(
                        f"{contig}, {position}, {nucleotide}: Singular matrix error, constant values, p-value set to 1; ",
                    )
                else:
                    handle_case(
                        nucleotide,
                        p_value=np.nan,
                        note=f"{nucleotide}: Singular matrix error occurred, p-value is NaN;",
                    )
                # Store warnings and continue
                position_result[f"{nucleotide}_warnings"] = (
                    "; ".join(warnings_list) or ""
                )
                position_result[f"{nucleotide}_n_warnings"] = len(warnings_list)
                # Skip remaining code for this freq_col
                continue

            # Dynamically extract the coefficient for `group`
            coef_name = next(
                (name for name in result.params.index if "group" in name), None
            )
            if coef_name is None:
                raise ValueError(
                    f"No group coefficient found for {freq_col}, {contig}, {position}"
                )

            # https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLMResults.html#statsmodels.regression.mixed_linear_model.MixedLMResults
            p_value = result.pvalues.get(coef_name)

            # Check for NaN p-values
            if pd.isna(p_value):
                # Post-hoc check: Only check perfect separation if we get NaN p-value
                if has_perfect_separation(sub_df, freq_col):
                    handle_case(
                        nucleotide,
                        note=f"{nucleotide}: perfect separation, p-value set to 0; ",
                        p_value=0,
                    )
                elif freq_values.nunique() == 1:
                    # Set p-value as 1 if there are no unique values for the nucleotide
                    handle_case(
                        nucleotide,
                        note=f"{nucleotide}: constant values, p-value set to 1; ",
                        p_value=1,
                    )
                else:
                    # Assign a p-value of NaN if none of the above two cases are true
                    var = freq_values.var()
                    handle_case(
                        nucleotide,
                        note=f"{nucleotide}: NaN p-value (variance: {var:.2e}); ",
                        p_value=np.nan,
                    )
            else:
                # Valid result case
                handle_case(
                    nucleotide,
                    note="",
                    p_value=p_value,
                    coef=result.params.get(coef_name),
                    t_value=result.tvalues.get(coef_name),
                )

            # Capture warnings regardless of outcome
            warnings_list = [str(warn.message) for warn in w]

        # Store warnings
        position_result[f"{nucleotide}_warnings"] = "; ".join(warnings_list) or ""
        position_result[f"{nucleotide}_n_warnings"] = len(warnings_list)

    return position_result


def main():
    logging.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Analyze allele frequency and perform significance tests.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_df",
        required=True,
        help="Path to input allele frequency dataframe (single) or mean changes dataframe (longitudnal).",
        type=str,
    )

    parser.add_argument(
        "--min_sample_num",
        type=int,
        default=4,
        help="Minimum number of samples per group required to perform the test.",
    )
    parser.add_argument(
        "--cpus",
        help="Number of processors to use.",
        type=int,
        default=cpu_count(),
    )
    parser.add_argument(
        "--outPath",
        help="Path to output file",
        type=str,
        required=True,
    )

    args = parser.parse_args()

    logging.info("Reading input dataframe...")
    df = pd.read_csv(args.input_df, sep="\t")

    grouped_positions = []
    for (contig, gene_id, position), sub_df in df.groupby(
        ["contig", "gene_id", "position"], dropna=False
    ):
        group_counts = sub_df.groupby("group").size()  # Count samples per group

        # Ensure both groups meet the min_sample_num threshold
        if (group_counts >= args.min_sample_num).sum() == 2:
            grouped_positions.append((contig, gene_id, position, sub_df))

    logging.info(
        f"Found {len(grouped_positions):,} unique (contig, gene_id, position) groups."
    )

    results = []
    logging.info(
        "Warnings will be captured and saved to the output table. warnings.simplefilter('default') is used."
    )

    # Use imap_unordered for parallel processing with a progress bar.
    with Pool(processes=args.cpus) as pool:
        for res in tqdm(
            pool.imap_unordered(run_model, grouped_positions),
            total=len(grouped_positions),
            desc="Running Linear Mixed Models",
        ):
            results.append(res)
    results_df = pd.DataFrame(results)

    initial_count = len(results_df)
    # Remove rows where all four nucleotide p-value columns are NaN
    results_df.dropna(
        subset=["A_p_value", "T_p_value", "G_p_value", "C_p_value"],
        how="all",
        inplace=True,
    )
    removed_count = initial_count - len(results_df)
    if removed_count:
        logging.info(f"Removed {removed_count:,} rows with all NaN p-values.")

    # Save the results to CSV.
    if args.outPath.endswith(".gz"):
        results_df.to_csv(args.outPath, index=False, sep="\t", compression="gzip")
    else:
        results_df.to_csv(args.outPath, index=False, sep="\t")
    logging.info(f"LMM modeling complete. Results saved to {args.outPath}")


if __name__ == "__main__":
    main()
