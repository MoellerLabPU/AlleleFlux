#!/usr/bin/env python3
"""Eligibility table: which genomes have enough replicates for which test.

QC writes one folder per timepoint with one ``{mag}_QC.tsv`` per genome, holding
per-sample facts only (did this sample pass the breadth and coverage
thresholds for this genome?).  This step reads the folders of ONE comparison
(two timepoints for longitudinal data, one for single) and, for the two tested
groups of that comparison, works out:

1. which subjects are readable at both timepoints (longitudinal only): only
   their samples count;
2. how many replicates each tested group has among those samples;
3. how many replicates the two tested groups share (the paired count);
4. which tests the genome is eligible for, given ``--min_sample_num``.

Small example (longitudinal, floor 2), groups A and B, replicate = cage:

    subject  group  cage  readable at pre  readable at post
    a1       A      c1    yes              yes
    a2       A      c2    yes              yes
    b1       B      c3    yes              yes
    b2       B      c4    yes              no

    -> a1, a2 and b1 count; b2 does not (not readable at post).
    -> A has 2 cages, B has 1; no cage is shared (paired count 0).
    -> unpaired test: not eligible (B is below the floor of 2).
"""

import argparse
import glob
import logging
import os
from typing import List, Optional

import numpy as np
import pandas as pd

from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.utilities import relabel_groups_from_metadata

logger = logging.getLogger(__name__)

# The label columns of a QC table this step needs.  They are read as TEXT:
# letting pandas guess turns a group named 20 into the number 20, and the
# group names given on the command line (always text) then match no row.
_QC_LABEL_COLUMNS = ("sample_id", "MAG_ID", "group", "subjectID", "replicate")
# The per-sample QC verdict: True when the sample passed the breadth AND the
# coverage threshold for this genome (the coverage check cascades from breadth).
_QC_VERDICT_COLUMN = "coverage_threshold_passed"
# QC folders a call must name: one per timepoint of the comparison.
_N_QC_DIRS = {"longitudinal": 2, "single": 1}
_QC_SUFFIX = "_QC.tsv"


def list_genomes(qc_dirs):
    """Every genome with a QC table in at least one of the folders, in name order.

    Name order (not the folder listing order, which depends on the file
    system) so that the eligibility table has the same row order on every run.

    Example: folders holding ``MAG_B_QC.tsv`` and ``MAG_A_QC.tsv`` ->
    ``["MAG_A", "MAG_B"]``.  The folder's summary files do not end in
    ``_QC.tsv`` and are not genomes.

    Raises ``FileNotFoundError`` when a folder does not exist.
    """
    genomes = set()
    for qc_dir in qc_dirs:
        if not os.path.isdir(qc_dir):
            raise FileNotFoundError(f"QC folder not found: {qc_dir}")
        for path in glob.glob(os.path.join(qc_dir, f"*{_QC_SUFFIX}")):
            genomes.add(os.path.basename(path)[: -len(_QC_SUFFIX)])
    return sorted(genomes)


def load_genome_qc(qc_dirs, mag_id):
    """One genome's per-sample QC rows from every folder, labels as text.

    Returns a frame with the columns ``sample_id, MAG_ID, group, subjectID,
    replicate, coverage_threshold_passed`` (True/False) and ``time`` when the
    tables have it; every other column of the QC tables is left unread.

    A header-only table (the genome has no profiled sample at that timepoint)
    contributes no row.  When every table is header-only the frame has the
    columns and no rows.

    Example: ``MAG_A`` has 4 rows in the ``pre`` folder and a header-only
    table in the ``post`` folder -> 4 rows, all with ``time == "pre"``.

    Raises:
        FileNotFoundError: a folder has no table for the genome.  Every genome
            has a table at every timepoint, so a missing one is a failed step.
        ValueError: a table lacks a needed column; the tables disagree on
            having ``time``; a label cell is blank; a sample is listed more
            than once; the verdict is not the text ``True`` or ``False``.
    """
    needed = set(_QC_LABEL_COLUMNS) | {_QC_VERDICT_COLUMN}
    header = None
    frames = []
    for qc_dir in qc_dirs:
        path = os.path.join(qc_dir, f"{mag_id}{_QC_SUFFIX}")
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"no QC table for {mag_id} in {qc_dir} (expected {path}): every "
                "genome has a table at every timepoint, header-only when it has "
                "no sample there, so a missing table means the QC step did not "
                "finish for this timepoint"
            )
        # Only the columns this step uses, all as text: a group named 20 or a
        # subject named 007 must stay the label it is in the sheet.
        df = pd.read_csv(
            path, sep="\t", dtype=str, usecols=lambda c: c in needed | {"time"}
        )
        missing = needed - set(df.columns)
        if missing:
            raise ValueError(f"{path}: QC table lacks column(s) {sorted(missing)}")
        if header is None:
            header = list(df.columns)
        elif set(df.columns) != set(header):
            raise ValueError(
                f"{path}: columns {sorted(df.columns)} differ from the other QC "
                f"table(s) of {mag_id} ({sorted(header)}); a time column must be "
                "in all of them or in none"
            )
        # A zero-row frame is left out of the concat: pandas would otherwise
        # widen the dtypes of the frames that do have rows.
        if not df.empty:
            frames.append(df[header])

    if not frames:
        empty = pd.DataFrame(columns=header)
        empty[_QC_VERDICT_COLUMN] = empty[_QC_VERDICT_COLUMN].astype(bool)
        return empty

    rows = pd.concat(frames, ignore_index=True)

    # A blank label would silently become a group or a subject of its own.
    blank = rows[header].isna().any(axis=1)
    if blank.any():
        raise ValueError(
            f"{mag_id}: blank cell in the QC table(s) for sample(s) "
            f"{sorted(rows.loc[blank, 'sample_id'].astype(str))}"
        )
    # One row per sample.  A sample listed twice means the same folder was
    # given twice or two folders overlap, and its subject would then look
    # like one seen at two timepoints.
    repeated = rows.loc[rows["sample_id"].duplicated(), "sample_id"]
    if not repeated.empty:
        raise ValueError(
            f"{mag_id}: sample(s) {sorted(set(repeated))} listed more than once "
            f"across the QC folders {list(qc_dirs)}"
        )
    # The verdict was read as text like everything else; turn it into real
    # True/False and refuse anything else instead of guessing.
    verdict = rows[_QC_VERDICT_COLUMN].map({"True": True, "False": False})
    if verdict.isna().any():
        raise ValueError(
            f"{mag_id}: {_QC_VERDICT_COLUMN} must be True or False, got "
            f"{sorted(set(rows.loc[verdict.isna(), _QC_VERDICT_COLUMN]))}"
        )
    rows[_QC_VERDICT_COLUMN] = verdict.astype(bool)
    return rows


def check_timepoints(df, data_type):
    """
    Check if MAGs (Metagenome-Assembled Genomes) are present in two timepoints for each subject.
    This function verifies that subjects have data from exactly two timepoints for longitudinal analysis.
    For single datatype, it skips the timepoint checks and sets two_timepoints_passed to match coverage_threshold_passed.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing sample data with 'subjectID', 'time', and 'coverage_threshold_passed' columns
    data_type : str
        Type of data analysis, either "single" or assumed to be longitudinal

    Returns
    -------
    pandas.DataFrame
        The input dataframe with an additional 'two_timepoints_passed' column indicating samples
        that come from subjects with exactly 2 timepoints

    Raises
    ------
    ValueError
        If more than 2 unique timepoints are found in the passed samples

    Notes
    -----
    - For longitudinal analysis, the function ensures there are exactly 2 unique timepoints globally
    - A sample is marked as passing if its subject has data at both timepoints and the sample has passed
      the coverage threshold
    - Logs information about subjects with valid timepoints
    """
    # For single datatype, skip timepoint checks
    if data_type == "single":
        # Set two_timepoints_passed to match coverage_threshold_passed for single data type
        df["two_timepoints_passed"] = df["coverage_threshold_passed"]
        return df

    # Initialize columns
    df["two_timepoints_passed"] = False

    # Get samples that passed coverage threshold
    passed_samples = df[df["coverage_threshold_passed"]]

    # Check that there are exactly 2 unique timepoints globally.
    unique_times = passed_samples["time"].unique()
    if len(unique_times) > 2:
        raise ValueError(f"More than 2 unique timepoints found: {unique_times}")

    if passed_samples.empty:
        logger.warning("No samples passed coverage threshold")
        return df

    # Find subjects with MAG in exactly 2 timepoints
    valid_subjects = (
        passed_samples.groupby("subjectID")["time"]
        .nunique()
        .where(lambda x: x == 2)
        .dropna()
        .index
    )

    # Mark passed samples from valid subjects
    valid_mask = passed_samples["subjectID"].isin(valid_subjects)
    df.loc[passed_samples[valid_mask].index, "two_timepoints_passed"] = True

    # Add debug info
    if not valid_subjects.empty:
        logger.debug(f"Subjects with 2 valid timepoints: {len(valid_subjects)}")
    else:
        logger.warning("No subjects have valid MAG in 2 timepoints")

    return df


def add_subject_count_per_group(df):
    """
    Add columns for subject count and replicate count per group to the DataFrame.

    This function calculates the number of unique subjects and replicates for each group
    in the dataset, but only considers samples that have passed timepoint validation
    (where 'two_timepoints_passed' is True). It then adds these counts as new columns
    to the original DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least the following columns:
        - 'two_timepoints_passed': Boolean indicating if the sample has passed timepoint validation
        - 'group': Group identifier
        - 'subjectID': Subject identifier
        - 'replicate': Replicate identifier

    Returns
    -------
    pandas.DataFrame
        The input DataFrame with two new columns added:
        - 'subjects_per_group': Number of unique subjects in the group (NaN for invalid samples)
        - 'replicates_per_group': Number of unique replicates in the group (NaN for invalid samples)

    Notes
    -----
    If no valid samples are found, both new columns will be set to NaN for all rows.
    The function logs a summary of subject and replicate counts per group for valid samples.
    """
    # Get samples that passed timepoint validation
    valid_subjects = df[df["two_timepoints_passed"]]

    if valid_subjects.empty:
        logger.debug("No valid subjects found")
        # Initialize the columns with NaN values
        df["subjects_per_group"] = np.nan
        df["replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "subjects_per_group"] = 0
        df.loc[df["two_timepoints_passed"], "replicates_per_group"] = 0
        return df

    # Count unique subjects per group in valid samples
    group_counts = (
        valid_subjects.groupby("group")["subjectID"]
        .nunique()  # Count unique subjects per group
        .reset_index(name="subjects_per_group")
    )

    # Count replicates per group in valid samples
    replicate_counts = (
        valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="replicates_per_group")
    )

    # Merge both counts into one counts DataFrame keyed by group.
    counts = group_counts.merge(replicate_counts, on="group", how="outer")

    # Merge the counts back into the original DataFrame based on group
    df = df.merge(counts, on="group", how="left")

    # For rows where two_timepoints_passed is False, set the counts to NaN.
    df.loc[
        ~df["two_timepoints_passed"], ["subjects_per_group", "replicates_per_group"]
    ] = np.nan

    logger.debug(f"Subject and replicate counts per group (valid samples):\n{counts}")
    return df


def count_paired_replicates(df):
    """
    Count the number of paired replicates per group in the dataset.

    This function identifies valid subjects (those that passed both timepoints) and counts
    replicates that appear in both groups. It ensures exactly two groups exist in the valid
    samples dataset.

    A "paired replicate" means the same replicate has samples in both experimental groups.
    This helps identify balanced experimental designs where the same experimental units
    are measured across different conditions.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing the following columns:
        - 'two_timepoints_passed': Boolean indicating if the sample has two valid timepoints
        - 'group': The experimental group identifier
        - 'replicate': The replicate identifier for the sample

    Returns
    -------
    pandas.DataFrame
        The original DataFrame with an additional column:
        - 'paired_replicates_per_group': Number of replicates per group that have samples in both groups.
          NaN for samples that failed timepoint check or whose replicate isn't paired.
          Zero (0) for samples that passed timepoint check but have no paired replicates.


    Notes
    -----
    - If there are no valid subjects at all, the function returns the DataFrame with
      'paired_replicates_per_group' column set to NaN for all rows.
    - If there are valid subjects but no paired replicates, the function sets
      'paired_replicates_per_group' to 0 for valid subjects (those that passed timepoint check).
    - For samples that failed the timepoint check, the column is always set to NaN.
    """
    valid_subjects = df[df["two_timepoints_passed"]]

    if valid_subjects.empty:
        logger.debug("No valid subjects found")
        # Initialize the column with NaN values
        df["paired_replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "paired_replicates_per_group"] = 0
        return df

    rep_group_count = (
        valid_subjects.groupby("replicate")["group"]
        .nunique()
        .reset_index(name="group_count")
    )
    paired_reps = rep_group_count[rep_group_count["group_count"] == 2]["replicate"]

    # Check if there are any paired replicates
    if paired_reps.empty:
        logger.warning("No paired replicates found across groups")
        # Initialize the column with NaN values
        df["paired_replicates_per_group"] = np.nan
        # Set 0 only for valid subjects that passed the timepoint check
        df.loc[df["two_timepoints_passed"], "paired_replicates_per_group"] = 0
        return df

    paired_valid_subjects = valid_subjects[
        valid_subjects["replicate"].isin(paired_reps)
    ]

    # Count the number of unique paired replicates per group
    paired_counts = (
        paired_valid_subjects.groupby("group")["replicate"]
        .nunique()
        .reset_index(name="paired_replicates_per_group")
    )

    # Merge the paired replicate counts back into the original DataFrame based on 'group'
    df = df.merge(paired_counts, on="group", how="left")

    # Set counts to NaN where:
    # - Samples failed timepoint check OR
    # - Samples passed timepoint check but replicate isn't paired
    valid_mask = df["two_timepoints_passed"]
    paired_mask = df["replicate"].isin(paired_reps)
    df["paired_replicates_per_group"] = np.where(
        valid_mask & paired_mask, df["paired_replicates_per_group"], np.nan
    )

    logger.debug(f"Replicates paired per group:\n{paired_counts}")
    return df


def pair_and_count(df, groups, data_type):
    """Mark the samples that count and count replicates over the tested groups.

    Three steps on one genome's per-sample rows (``load_genome_qc``):

    1. ``two_timepoints_passed``: longitudinal, a sample counts when it passed
       QC and its subject has a passing sample at BOTH timepoints; single, a
       sample counts when it passed QC.  Decided over every group's rows, as a
       subject is a subject whatever its group.
    2. Keep the rows of the two tested groups (``groups``; None keeps all).
    3. ``subjects_per_group`` / ``replicates_per_group``: unique subjects and
       replicates of each group among its counting samples.
       ``paired_replicates_per_group``: replicates that both tested groups
       have among their counting samples.

    Step 3 comes after step 2 on purpose: "shared by the tested pair" can only
    be asked of the two tested groups.  Counted over every group of a run,
    three groups that reuse the same replicate IDs have no replicate "in
    exactly two groups", and the paired test was never eligible.

    Example (longitudinal, tested pair A-B; block r1 has one subject in each
    of A, B and C, all readable at both timepoints): after step 2 the C rows
    are gone; A and B each have 1 replicate and 1 paired replicate.

    Returns the tested groups' rows with the four columns added; counts are
    NaN on rows that do not count.

    Raises ``ValueError`` when longitudinal rows have no ``time`` column or
    span more than two timepoints.
    """
    if data_type == "longitudinal" and "time" not in df.columns:
        raise ValueError(
            "longitudinal eligibility needs a time column in the QC tables; "
            f"columns found: {list(df.columns)}"
        )
    df = check_timepoints(df, data_type)
    if groups is not None:
        # .copy(): the counting functions add columns, and a filtered view of
        # the full frame would raise pandas' chained-assignment warning.
        df = df[df["group"].isin(groups)].copy()
    df = add_subject_count_per_group(df)
    df = count_paired_replicates(df)
    return df


def process_data(df, unique_groups, mag_id, min_sample_num, data_type):
    """
    Process and summarize data eligibility based on replicate counts for each group.

    This function checks for consistency in replicate and paired replicate counts within each given group.
    It computes a summary pivot table containing the number of replicates per group and determines eligibility
    for unpaired and paired tests as well as single sample eligibility (dependent on the data type).

    Parameters:
        df (pandas.DataFrame): DataFrame containing the data with at least the columns "group",
                               "replicates_per_group", and "paired_replicates_per_group".
        unique_groups (list): List of unique group identifiers present in the DataFrame.
        mag_id (Any): Identifier of the current MAG (or sample) being processed.
        min_sample_num (int): The minimum number of replicates required for eligibility.
        data_type (str): Type of the data; if set to "longitudinal", single sample eligibility is assessed
                         per group. Otherwise, the single sample eligibility values are set to NaN.

    Returns:
        dict: A dictionary containing the pivoted summary of replicates and eligibility flags. Keys include:
              - 'unpaired_test_eligible': A boolean indicating if both groups have replicates greater than or equal to min_sample_num.
              - 'paired_test_eligible': A boolean indicating if both groups have paired replicates greater than or equal to min_sample_num.
              - 'single_sample_eligible_<group>': For each group, a boolean (if data_type == "longitudinal") showing
                if the replicates per group are greater than or equal to min_sample_num, or NaN otherwise.

    Raises:
        ValueError: If inconsistent replicate counts are detected within any group, indicating an error for the given MAG.
    """
    # Verify that replicates_per_group and paired_replicates_per_group are consistent within each group
    for group in unique_groups:
        group_df = df[df["group"] == group]
        if (
            group_df["replicates_per_group"].nunique() != 1
            or group_df["paired_replicates_per_group"].nunique() != 1
        ):
            raise ValueError(
                f"MAG {mag_id}, group {group}: Inconsistent replicates detected. Skipping."
            )

    summary = df.groupby(["MAG_ID", "group"], as_index=False).agg(
        replicates_per_group=("replicates_per_group", "first"),
        paired_replicates_per_group=("paired_replicates_per_group", "first"),
    )

    # Pivot the table so that each row is a unique combination of MAG_ID and the groups become separate columns.
    pivot = summary.pivot(index=["MAG_ID"], columns="group")
    pivot.columns = [f"{stat}_{group}" for stat, group in pivot.columns]
    pivot = pivot.reset_index()

    g1, g2 = unique_groups

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
        for group in unique_groups:
            pivot[f"single_sample_eligible_{group}"] = (
                pivot[f"replicates_per_group_{group}"] >= min_sample_num
            )
    elif data_type == "single":
        # Set single sample eligibility columns to NaN for single timepoint data
        for group in unique_groups:
            pivot[f"single_sample_eligible_{group}"] = np.nan

    return pivot.to_dict(orient="records")[0]


def analyze_qc_dirs(
    qc_dirs: List[str],
    min_sample_num: int,
    data_type: str,
    groups: Optional[List[str]] = None,
    permuted_metadata: Optional[str] = None,
) -> pd.DataFrame:
    """Build the eligibility table of one comparison from its per-timepoint QC folders.

    For every genome with a QC table: read its rows from each folder
    (``load_genome_qc``), relabel the groups for a permuted run, mark the
    samples that count and count replicates over the tested groups
    (``pair_and_count``), then decide the tests (``process_data``).

    A genome reaches the table only when BOTH tested groups have at least one
    counting sample; the others are logged and left out.

    Parameters:
        qc_dirs: the QC folders, one per timepoint of the comparison: exactly
            two for ``longitudinal``, exactly one for ``single``.
        min_sample_num: the replicate floor of every test.
        data_type: ``"longitudinal"`` or ``"single"``.
        groups: the two tested groups.  None takes the groups found in each
            genome's rows (the genome is left out unless there are exactly two).
        permuted_metadata: optional permuted sheet of a null run; the groups
            are relabelled from it (joined on ``sample_id``) BEFORE anything
            is counted, so the counts are those of the permuted labels.

    Returns:
        One row per eligible-or-not genome, in genome name order, with the
        columns ``MAG_ID``, ``replicates_per_group_<g>`` and
        ``paired_replicates_per_group_<g>`` for both groups,
        ``unpaired_test_eligible``, ``paired_test_eligible`` and
        ``single_sample_eligible_<g>`` (the tested pair in the order given).

    Example: two genomes, tested pair A-B, floor 4.  ``MAG_A`` has 4 counting
    cages in each group and none shared -> replicates 4 and 4, paired 0 and 0,
    unpaired eligible, paired not.  ``MAG_B`` has samples at one timepoint
    only -> no counting sample, left out.

    Raises:
        ValueError: ``groups`` is not two names; the number of folders does
            not fit the data type; no QC table in the folders; no genome
            reaches the table.
        FileNotFoundError: a folder is missing, or lacks a genome's table.
    """
    if groups is not None and len(groups) != 2:
        raise ValueError(
            f"--groups must specify exactly 2 group names for eligibility filtering, "
            f"got {len(groups)}: {groups}"
        )
    if data_type not in _N_QC_DIRS:
        raise ValueError(
            f"data_type must be one of {sorted(_N_QC_DIRS)}, got {data_type!r}"
        )
    expected = _N_QC_DIRS[data_type]
    if len(qc_dirs) != expected:
        raise ValueError(
            f"{data_type} data need exactly {expected} QC folder(s), one per "
            f"timepoint of the comparison; got {len(qc_dirs)}: {list(qc_dirs)}"
        )

    genomes = list_genomes(qc_dirs)
    if not genomes:
        raise ValueError(f"No QC files found in {list(qc_dirs)}")

    results = []
    for mag_id in genomes:
        df = load_genome_qc(qc_dirs, mag_id)
        if df.empty:
            logger.info(
                f"MAG {mag_id}: no sample at the timepoint(s) of this comparison. "
                "Skipping."
            )
            continue

        # Permuted (null) run: put the permuted group labels on the rows before
        # anything is counted.  Only `group` moves; subjects, replicates and
        # timepoints stay with their samples.
        if permuted_metadata:
            df = relabel_groups_from_metadata(df, permuted_metadata)

        df = pair_and_count(df, groups, data_type)
        passed_df = df[df["two_timepoints_passed"]]

        # The tested groups that have at least one counting sample, in a fixed
        # order: the pair as given, or name order when no pair was given.  (The
        # order the groups happen to appear in the QC rows depends on which
        # sample a QC worker finished first, so it must not decide the order of
        # the table's columns.)
        present = set(passed_df["group"])
        if groups is not None:
            tested = [g for g in groups if g in present]
        else:
            tested = sorted(present)

        if len(tested) != 2:
            logger.warning(
                f"MAG {mag_id}: Found {len(tested)} groups after filtering "
                f"(expected exactly 2). Skipping."
            )
            continue

        logger.info(f"Processing MAG: {mag_id}")
        results.append(
            process_data(passed_df, tested, mag_id, min_sample_num, data_type)
        )

    if not results:
        raise ValueError("No MAGs passed eligibility criteria.")

    return pd.DataFrame(results)


def main():
    setup_logging()

    parser = argparse.ArgumentParser(
        description="Generate test eligibility summary table from QC files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--qc_dirs",
        required=True,
        nargs="+",
        metavar="DIR",
        help="QC folders holding the *_QC.tsv files, one folder per timepoint of "
        "the comparison: exactly two for longitudinal data, exactly one for single.",
    )
    parser.add_argument(
        "--min_sample_num", type=int, default=4, help="Minimum replicates required"
    )
    parser.add_argument(
        "--output_file", required=True, help="Output file for summary table"
    )
    parser.add_argument(
        "--data_type",
        choices=["single", "longitudinal"],
        default="longitudinal",
        help="Type of data: single timepoint or longitudinal",
    )
    parser.add_argument(
        "--groups",
        nargs=2,
        default=None,
        metavar=("GROUP1", "GROUP2"),
        help="The two tested groups. The QC folders cover every group of the run; "
        "replicates and paired replicates are counted over these two only.",
    )
    parser.add_argument(
        "--permuted_metadata",
        type=str,
        default=None,
        help=(
            "Optional path to a permuted metadata TSV (null/control run).  When "
            "given, each sample's group is re-derived from this sheet (joined on "
            "sample_id) before the replicates are counted; no profile reads.  "
            "Omit for a normal run."
        ),
    )

    args = parser.parse_args()

    eligibility_table = analyze_qc_dirs(
        args.qc_dirs,
        args.min_sample_num,
        args.data_type,
        groups=args.groups,
        permuted_metadata=args.permuted_metadata,
    )
    eligibility_table.to_csv(args.output_file, sep="\t", index=False)

    logger.info(f"Eligibility summary table written to {args.output_file}")


if __name__ == "__main__":
    main()
