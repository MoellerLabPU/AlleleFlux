#!/usr/bin/env python
import argparse
import logging
import warnings
from multiprocessing import Pool, cpu_count

import pandas as pd
import statsmodels.formula.api as smf

# import supress_warning
from tqdm import tqdm


def run_model(args):
    freq_cols = ["A_frequency", "T_frequency", "G_frequency", "C_frequency"]

    contig, position, sub_df = args

    position_result = {
        "contig": contig,
        "position": position,
        "n_samples": len(sub_df),
    }

    for freq_col in freq_cols:
        nucleotide = freq_col.split("_")[0]
        warnings_list = []
        with warnings.catch_warnings(record=True) as w:
            # https://docs.python.org/3/library/warnings.html#the-warnings-filter
            warnings.simplefilter("default")
            # Fit Linear Mixed Model
            model = smf.mixedlm(
                f"{freq_col} ~ group", data=sub_df, groups=sub_df["replicate"]
            )
            result = model.fit()
            # Dynamically extract the coefficient for `group`
            coef_names = [name for name in result.params.index if "group" in name]
            if coef_names:
                coef_name = coef_names[0]
                # https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLMResults.html#statsmodels.regression.mixed_linear_model.MixedLMResults
                position_result[f"{nucleotide}_coef"] = result.params.get(
                    coef_name, None
                )
                position_result[f"{nucleotide}_p_value"] = result.pvalues.get(
                    coef_name, None
                )
                # Get the t-value for the coefficient
                position_result[f"{nucleotide}_t_value (z-score)"] = result.tvalues.get(
                    coef_name, None
                )
            else:
                raise ValueError(f"Could not find coefficient for {freq_col}")
            # Append any captured warnings to warnings_list.
            for warn in w:
                warnings_list.append(str(warn.message))
        position_result[f"{nucleotide}_warnings"] = (
            "; ".join(warnings_list) if warnings_list else ""
        )
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
    for (contig, position), sub_df in df.groupby(["contig", "position"]):
        group_counts = sub_df.groupby("group").size()  # Count samples per group

        # Ensure both groups meet the min_sample_num threshold
        if (
            group_counts
            >= args.min_sample_num
            # Ensure both groups meet criteria
        ).sum() == 2:
            grouped_positions.append((contig, position, sub_df))

    logging.info(f"Found {len(grouped_positions):,} unique (contig, position) groups.")

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

    # Save the results to CSV.
    results_df.to_csv(args.outPath, index=False)
    logging.info(f"LME modeling complete. Results saved to {args.outPath}")


if __name__ == "__main__":
    main()
