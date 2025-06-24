#!/usr/bin/env python
import argparse
import logging
import os
import re
import warnings
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from tqdm import tqdm

from alleleflux.scripts.utilities.utilities import load_and_filter_data

# Accept any suffix that does *not* begin with "diff"
_FREQ_PATTERN = re.compile(r"^[ATGC]_frequency_(?!diff\b).*")

NUCLEOTIDES = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]


def plan_csv_loading(dtype_map, input_time_format, df_path):
    """
    Plan for loading CSV files with specific data types.
    This function is used to define the expected data types for each column
    in the input DataFrame.
    """
    if input_time_format == "column":
        dtype_map["time"] = str
        dtype_map.update({nuc: "float32" for nuc in NUCLEOTIDES})

    header = pd.read_csv(df_path, sep="\t", nrows=0).columns

    if input_time_format == "suffix":
        for col in header:
            if _FREQ_PATTERN.match(col):
                dtype_map[col] = "float32"

    return dtype_map


def convert_suffix_to_long(df_source):
    """
    Convert wide-format frequency columns with timepoint suffixes into long format.
    Expects columns like A_frequency_<suffix>, T_frequency_<suffix>, etc.
    Returns a DataFrame with a new 'time' column and base stubname columns.

    Parameters
    ----------
    df_source : str | pandas.DataFrame
        Either a path to a TSV file OR a DataFrame already in memory.

    Returns
    -------
    pandas.DataFrame

    """
    # Input is a string
    if isinstance(df_source, str):
        df = pd.read_csv(df_source, sep="\t")
    else:
        # Input is already a DataFrame
        df = df_source.copy()

    # Identify frequency columns with timepoint suffixes
    suffixed_cols = [c for c in df.columns if _FREQ_PATTERN.match(c)]
    if not suffixed_cols:
        raise ValueError(
            "No suffixed frequency columns detected for suffix-based time format."
        )
    # Extract unique suffixes (timepoints)
    suffixes = sorted({c.split("_frequency_", 1)[1] for c in suffixed_cols})
    if len(suffixes) < 2:
        raise ValueError(
            f"Less than two timepoint suffixes found in columns: {suffixes}"
        )

    # Define id_vars for pivoting
    id_vars = [c for c in df.columns if c not in suffixed_cols]

    # Perform wide-to-long conversion
    df_long = pd.wide_to_long(
        df, stubnames=NUCLEOTIDES, i=id_vars, j="time", sep="_", suffix=r".+"
    ).reset_index()

    logging.info(
        f"Converted wide-format frequencies to long format with timepoints: {suffixes}"
    )
    return df_long


def has_perfect_separation(sub_df, freq_col, predictor_col):
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
        The sub_df is analyzed based on groupings from predictor_col

    freq_col : str
        Name of the column to analyze for perfect separation within groups.

    Returns
    -------
    bool
        True if perfect separation is detected, False otherwise.
    """

    # Group by predictor_col and compute mean and variance
    group_stats = sub_df.groupby(predictor_col)[freq_col].agg(["mean", "var"])

    # Check:
    # 1. All groups have zero variance (var == 0)
    # 2. At least two groups have different means
    return (group_stats["var"] == 0).all() and (group_stats["mean"].nunique() > 1)


def detect_column_suffix(data_type):
    # Identify if the columns are mean changes or raw frequencies
    base_freq_cols = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]
    if data_type == "longitudinal":
        freq_cols = [f"{nuc}_diff_mean" for nuc in base_freq_cols]
    elif data_type == "single" or data_type == "across_time":
        freq_cols = base_freq_cols
    else:
        raise ValueError(f"Unknown data_type: {data_type}")
    return freq_cols


def run_model(args):
    contig, gene_id, position, sub_df, data_type = args

    position_result = {
        "contig": contig,
        "gene_id": gene_id,
        "position": position,
        "num_pairs": len(sub_df),
        "notes": "",
    }

    predictor_col = "group"
    if data_type == "across_time":
        predictor_col = "time"
        # If analyzing across time for a specific group, store that group.
        if "group" in sub_df.columns and sub_df["group"].nunique() == 1:
            position_result["group_analyzed"] = sub_df["group"].iloc[0]

    def handle_case(nucleotide, note, p_value=np.nan, coef=np.nan, t_value=np.nan):
        position_result[f"{nucleotide}_coef_LMM"] = coef
        position_result[f"{nucleotide}_p_value_LMM"] = p_value
        position_result[f"{nucleotide}_t_value_LMM (z-score)"] = t_value
        position_result["notes"] += note

    # Detect the suffix of the columns based on the data type
    # This is used to determine if the columns are mean changes or raw frequencies
    freq_cols = detect_column_suffix(data_type)

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
                    f"{freq_col} ~ {predictor_col}",
                    sub_df,
                    groups=sub_df["replicate"],  # Use dynamic predictor_col
                )
                result = model.fit()

            except np.linalg.LinAlgError as e:
                # Capture warnings first before continuing
                warnings_list = [str(warn.message) for warn in w]
                # These checks are done after model fitting to capture the warnings
                if has_perfect_separation(
                    sub_df, freq_col, predictor_col
                ):  # Pass predictor_col
                    handle_case(
                        nucleotide,
                        note=f"{nucleotide}: Singular matrix error, perfect separation with {predictor_col}, p-value set to 0; ",
                        p_value=0,
                    )
                    logging.info(
                        f"{contig}, {position}, {nucleotide}: Singular matrix error, perfect separation with {predictor_col}, p-value set to 0; ",
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

            # Dynamically extract the coefficient for the predictor_col
            # The coefficient name might be like 'group[T.group_name]' or 'time[T.timepoint_name]'
            coef_name = next(
                (name for name in result.params.index if predictor_col in name), None
            )
            if coef_name is None:
                raise ValueError(
                    f"No {predictor_col} coefficient found for {freq_col}, {contig}, {position}. Params: {result.params.index}"
                )

            # https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLMResults.html#statsmodels.regression.mixed_linear_model.MixedLMResults
            p_value = result.pvalues.get(coef_name)

            # Check for NaN p-values
            if pd.isna(p_value):
                # Post-hoc check: Only check perfect separation if we get NaN p-value
                if has_perfect_separation(
                    sub_df, freq_col, predictor_col
                ):  # Pass predictor_col
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
        format="[%(asctime)s %(levelname)s] %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logging.INFO,
    )

    parser = argparse.ArgumentParser(
        description="Analyze allele frequency and perform significance tests.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input_df",
        required=True,
        help="Path to input allele frequency dataframe (single, across_time) or mean changes dataframe (longitudnal).",
        type=str,
    )
    parser.add_argument(
        "--preprocessed_df",
        help="Optional path to a preprocessed dataframe to filter positions. Only used when --data_type is 'across_time'.",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--data_type",
        help="Type of data to analyze: longitudinal, single, or across_time",
        type=str,
        choices=["longitudinal", "single", "across_time"],  # Added across_time
        default="longitudinal",
    )
    parser.add_argument(
        "--group_to_analyze",
        help="For 'across_time' data_type only: name of the group to analyze across timepoints. This group will be fixed, and comparison will be between timepoints.",
        type=str,
        default=None,  # Default to None
    )
    parser.add_argument(
        "--input_time_format",
        help="Format of time information in the input dataframe. 'column': expects a 'time' column. 'suffix': expects timepoints as suffixes in frequency columns (e.g., A_frequency_t1, C_frequency_t2). 'suffix' is only valid for 'across_time' data_type, for the rest 'column' is used.",
        type=str,
        choices=["column", "suffix"],
        default="column",
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
        "--output_dir",
        help="Path to output directory.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--mag_id",
        help="MAG ID to process",
        type=str,
        required=True,
    )

    args = parser.parse_args()
    # Validate input arguments
    if args.input_time_format == "suffix" and args.data_type != "across_time":
        parser.error(
            "--input_time_format 'suffix' is only valid with --data_type 'across_time'."
        )

    logging.info("Reading input dataframe...")
    dtype_map = {
        "subjectID": str,
        "contig": str,
        "gene_id": str,
        "position": int,
        "group": str,  # Use category for group to save memory
        "replicate": str,  # Use category for replicate to save memory
    }
    if args.preprocessed_df and args.data_type == "across_time":
        # If preprocessed_df is provided, load and filter input_df based on it
        logging.info(
            f"Preprocessed dataframe provided. Loading and filtering input_df based on {args.preprocessed_df}"
        )

        df = load_and_filter_data(
            input_df_path=args.input_df,
            preprocessed_df_path=args.preprocessed_df,
            mag_id=args.mag_id,
            dtype_map=plan_csv_loading(
                dtype_map, args.input_time_format, args.input_df
            ),
            group_to_analyze=(args.group_to_analyze),
        )
        if df.empty:
            raise ValueError(
                f"DataFrame is empty after filtering with preprocessed_df for MAG {args.mag_id}. No LMM analysis will be performed."
            )

        if args.input_time_format == "suffix":
            logging.info(
                "Converting frequency columns with timepoint suffixes to long format."
            )
            df = convert_suffix_to_long(df)
    else:
        if args.data_type == "across_time":
            logging.info(
                f"No preprocessed dataframe provided. Loading input_df directly from {args.input_df}"
            )
            if args.input_time_format == "suffix":
                logging.info(
                    "Converting frequency columns with timepoint suffixes to long format."
                )
                df = convert_suffix_to_long(args.input_df)
            else:
                dtype_map = plan_csv_loading(
                    dtype_map, args.input_time_format, args.input_df
                )

                df = pd.read_csv(
                    args.input_df,
                    sep="\t",
                    dtype=dtype_map,
                    usecols=dtype_map.keys(),
                    memory_map=True,
                )
        else:
            df = pd.read_csv(
                args.input_df,
                sep="\t",
            )
        logging.info(f"Loaded {df.shape[0]:,} rows from {args.input_df}.")

    # Initial filtering for across_time mode
    if args.data_type == "across_time":
        if not args.group_to_analyze:
            parser.error("--group_to_analyze is required for --data_type across_time")

        if args.group_to_analyze not in df["group"].unique():
            raise ValueError(
                f"Group '{args.group_to_analyze}' not found in input data. LMM analysis will yield no results."
            )
        df = df[df["group"] == args.group_to_analyze].copy()
        logging.info(f"Filtered data for group: {args.group_to_analyze}")

        # Handle time format based on input_time_format
        if "time" not in df.columns:
            raise ValueError("'time' column not found in input data. Exiting.")

        # Ensure at least two timepoints for analysis
        if df["time"].nunique() < 2:
            raise ValueError(
                f"Less than two unique timepoints found for group '{args.group_to_analyze}'. Unique timepoints: {df['time'].unique()}"
            )

    grouped_positions = []
    for (contig, gene_id, position), sub_df_orig in df.groupby(
        ["contig", "gene_id", "position"], dropna=False
    ):
        sub_df = sub_df_orig.copy()  # Work with a copy to avoid SettingWithCopyWarning

        if args.data_type == "across_time":
            unique_times = sub_df["time"].unique()
            if len(unique_times) != 2:
                # logging.debug(f"Skipping {contig}-{position} for group {args.group_to_analyze}: Expected 2 timepoints, found {len(unique_times)} ({unique_times})")
                continue

            # Count unique replicates per timepoint for the current contig/position/group
            time_point_replicate_counts = sub_df.groupby("time")["replicate"].nunique()

            if not (
                (time_point_replicate_counts >= args.min_sample_num).all()
                and len(time_point_replicate_counts) == 2
            ):
                # logging.debug(f"Skipping {contig}-{position} for group {args.group_to_analyze}: Not enough unique replicates per timepoint. Counts: {time_point_replicate_counts.to_dict()}, Min required: {args.min_sample_num}")
                continue
        else:  # Original logic for 'single' and 'longitudinal'
            uniques = sub_df["group"].unique()
            if len(uniques) != 2:
                # logging.debug(f"Skipping {contig}-{position}: Expected 2 groups, found {len(uniques)} ({uniques})")
                continue

            group_replicate_counts = sub_df.groupby("group")["replicate"].nunique()
            if not (
                (group_replicate_counts >= args.min_sample_num).all()
                and len(group_replicate_counts) == 2
            ):
                # logging.debug(f"Skipping {contig}-{position}: Not enough unique replicates per group. Counts: {group_replicate_counts.to_dict()}, Min required: {args.min_sample_num}")
                continue

        grouped_positions.append((contig, gene_id, position, sub_df, args.data_type))

    logging.info(
        f"Found {len(grouped_positions):,} unique (contig, gene_id, position) groups with sufficient replicates."
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
        subset=["A_p_value_LMM", "T_p_value_LMM", "G_p_value_LMM", "C_p_value_LMM"],
        how="all",
        inplace=True,
    )
    removed_count = initial_count - len(results_df)
    if removed_count:
        logging.info(f"Removed {removed_count:,} rows with all NaN p-values.")

    output_file_name = f"{args.mag_id}_lmm.tsv.gz"
    if args.data_type == "across_time" and args.group_to_analyze:
        output_file_name = (
            f"{args.mag_id}_lmm_across_time_{args.group_to_analyze}.tsv.gz"
        )

    output_path = os.path.join(args.output_dir, output_file_name)

    os.makedirs(args.output_dir, exist_ok=True)
    logging.info(f"Saving LMM results for MAG {args.mag_id} to {output_path}")
    results_df.to_csv(
        output_path,
        index=False,
        sep="\t",
        compression="gzip",
    )

    logging.info(f"LMM modeling complete. Results saved to {output_path}")


if __name__ == "__main__":
    main()
