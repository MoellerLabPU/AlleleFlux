#!/usr/bin/env python3

import argparse
import re

# Regex to capture groups, mag, timepoints
# Format:
#   significance_analysis_per_mag-groups=...,mag=...,timepoints=...
# significance_analysis_per_mag-groups=fat_control,mag=SLG1009_DASTool_bins_111,timepoints=end_post-61553487.out
# significance_analysis_per_mag-groups=fat_control,mag=SLG1122_DASTool_bins_78,timepoints=end_post-61553782.out
# significance_analysis_per_mag-groups=fat_control,mag=SLG992_DASTool_bins_103,timepoints=end_post-61552666.out
PATTERN = re.compile(
    r"^significance_analysis_per_mag-groups=([^,]+),mag=([^,]+),timepoints=(.+)$"
)


def parse_line(line):
    line = line.strip()
    match = PATTERN.match(line)
    if match:
        groups, mag, timepoints = match.groups()
        return groups.strip(), mag.strip(), timepoints.strip()
    return None


def make_table(input_file, output_file):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        # Print header
        fout.write("mag\tgroups\ttimepoints\n")

        for line in fin:
            parsed = parse_line(line)
            groups, mag, timepoints = parsed

            # Split timepoints at the first dash and keep only the first part
            # e.g. "end_post-61553487.out" => "end_post"
            # or "pre_end-61552644.out" => "pre_end"
            timepoints_short = timepoints.split("-", 1)[0].strip()

            # Write in the order: mag, groups, timepoints
            fout.write(f"{mag}\t{groups}\t{timepoints_short}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert lines from significance_analysis_per_mag format into a blacklist TSV with mag, groups, timepoints."
    )
    parser.add_argument(
        "--input_file",
        help="Text file with lines: significance_analysis_per_mag-groups=...,mag=...,timepoints=...",
        required=True,
    )
    parser.add_argument(
        "--output_file", help="Path to the output .tsv file.", required=True
    )
    args = parser.parse_args()
    make_table(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
