#!/usr/bin/env python3
"""
Terminal nucleotide analysis for specified group and timepoint samples.

This script processes genomic data to identify terminal nucleotides in specified samples
by analyzing nucleotide frequency profiles for significant sites.

Three methods for determining terminal nucleotide:
1. Mean frequency: Nucleotide with highest mean frequency at terminal timepoint
2. Majority voting: Most common major allele across terminal samples
3. Frequency change (optional): Nucleotide with greatest positive frequency change
   from initial to terminal timepoint, calculated per-subject then averaged

For frequency change analysis, samples are paired by subjectID (same individual at
initial vs terminal timepoint) to compute within-subject changes before averaging.
"""

import argparse
import logging
import math
import multiprocessing
from collections import Counter
from functools import partial
from pathlib import Path
from typing import Optional

import pandas as pd
from tqdm import tqdm

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# Global constant for nucleotide symbols
NUCLEOTIDES = ["A", "C", "G", "T"]

# Standard dtypes for profile files
PROFILE_DTYPES = {
    "contig": str,
    "position": int,
    "total_coverage": float,
    **{nuc: "int32" for nuc in NUCLEOTIDES},
}


def find_profile_files(profile_dir, sample_ids, mag_ids):
    """
    Find profile files in subdirectories by sample ID and MAG ID.

    Parameters:
    -----------
    profile_dir : str
        Root directory containing profile subdirectories by sample ID
    sample_ids : list
        List of sample IDs to look for
    mag_ids : list
        List of MAG IDs to filter profile files by. If None, uses wildcard pattern.

    Returns:
    --------
    dict : Dictionary mapping sample_id to list of profile file paths
    """
    profile_files = {}
    profile_dir_path = Path(profile_dir)

    if not profile_dir_path.exists():
        raise FileNotFoundError(f"Profile directory not found: {profile_dir}")

    logger.info(
        f"Searching for profile files for {len(sample_ids)} sample(s) and {len(mag_ids)} MAG(s)"
    )

    # Check if profile_dir contains subdirectories (new structure) or files directly
    subdirs = [d for d in profile_dir_path.iterdir() if d.is_dir()]

    if subdirs:
        # New structure: subdirectories by sample_id
        logger.info(f"Found {len(subdirs)} subdirectories in profile directory")
        for subdir in subdirs:
            sample_id = subdir.name
            if sample_id in sample_ids:
                for mag_id in mag_ids:
                    base = f"{sample_id}_{mag_id}_profiled.tsv"
                    candidates = list(subdir.glob(base)) + list(
                        subdir.glob(base + ".gz")
                    )
                    if candidates:
                        profile_files.setdefault(sample_id, []).extend(
                            map(str, candidates)
                        )
    # Log summary
    total_files_found = sum(len(files) for files in profile_files.values())
    samples_with_files = len(profile_files)

    logger.info(
        f"Profile file search complete: found {total_files_found} files for {samples_with_files}/{len(sample_ids)} samples"
    )

    return profile_files


def get_terminal_samples(metadata_file, group, timepoint):
    """
    Get sample IDs that match the specified group and timepoint.

    Parameters:
    -----------
    metadata_file : str
        Path to metadata file
    group : str
        Group name to filter for
    timepoint : str
        Timepoint to filter for

    Returns:
    --------
    list : List of sample IDs matching the criteria
    """
    logger.info(
        f"Loading metadata to identify samples for group '{group}' and timepoint '{timepoint}'"
    )
    metadata = pd.read_csv(metadata_file, sep="\t")

    if (
        "group" not in metadata.columns
        or "time" not in metadata.columns
        or "sample_id" not in metadata.columns
    ):
        raise ValueError(
            "Metadata must contain 'group', 'time', and 'sample_id' columns"
        )

    # Filter for matching group and timepoint
    matching_samples = metadata[
        (metadata["group"] == group) & (metadata["time"] == timepoint)
    ]["sample_id"].tolist()

    if not matching_samples:
        raise ValueError(
            f"No samples found for group '{group}' and timepoint '{timepoint}' in metadata"
        )

    logger.info(
        f"Found {len(matching_samples)} matching samples: {matching_samples[:5]}{'...' if len(matching_samples) > 5 else ''}"
    )

    return matching_samples


def get_subject_sample_mapping(
    metadata_file, group, terminal_timepoint, initial_timepoint
):
    """
    Create mapping of subjects to their sample IDs at each timepoint.

    Parameters:
    -----------
    metadata_file : str
        Path to metadata file
    group : str
        Group name to filter for
    terminal_timepoint : str
        Terminal timepoint identifier
    initial_timepoint : str
        Initial timepoint identifier

    Returns:
    --------
    dict : Dictionary mapping subject_id to dict with 'terminal' and 'initial' sample_ids
    """
    metadata = pd.read_csv(metadata_file, sep="\t")

    required_cols = ["group", "time", "sample_id", "subjectID"]
    missing = set(required_cols) - set(metadata.columns)
    if missing:
        raise ValueError(
            f"Metadata must contain {required_cols} columns for frequency change analysis. "
            f"Missing: {missing}"
        )

    # Filter for the specified group
    group_metadata = metadata[metadata["group"] == group]

    # Build subject-to-sample mapping
    subject_mapping = {}

    for _, row in group_metadata.iterrows():
        subject = row["subjectID"]
        sample_id = row["sample_id"]
        timepoint = row["time"]

        if subject not in subject_mapping:
            subject_mapping[subject] = {"terminal": None, "initial": None}

        if timepoint == terminal_timepoint:
            subject_mapping[subject]["terminal"] = sample_id
        elif timepoint == initial_timepoint:
            subject_mapping[subject]["initial"] = sample_id

    # Filter to only subjects with both timepoints
    complete_subjects = {
        subj: mapping
        for subj, mapping in subject_mapping.items()
        if mapping["terminal"] is not None and mapping["initial"] is not None
    }

    if not complete_subjects:
        raise ValueError(
            f"No subjects found with samples at both timepoints "
            f"('{initial_timepoint}' and '{terminal_timepoint}') in group '{group}'"
        )

    logger.info(
        f"Found {len(complete_subjects)} subjectID(s) with paired samples for frequency change analysis"
    )

    return complete_subjects


def perform_majority_voting(major_alleles):
    """
    Perform majority voting to determine terminal nucleotide from major alleles.

    When a sample has tied major alleles, all tied alleles contribute fractionally
    to the vote count (each gets 1/n_tied votes) to avoid biasing toward any single nucleotide.

    Parameters:
    -----------
    major_alleles : list
        List of major allele strings from all samples, or None values for missing data.
        Each string can be a single nucleotide or comma-separated tied nucleotides.

    Returns:
    --------
    str : The terminal nucleotide(s) determined by majority voting (comma-separated if tied)
    """
    # Count votes with fractional contribution for ties
    allele_counts = Counter()

    for major_allele in major_alleles:
        if major_allele is None:
            continue
        # major_allele can be single or comma-separated nucleotides
        votes = major_allele.split(",")
        # Each tied nucleotide gets equal fractional vote
        vote_weight = 1.0 / len(votes)

        for nucleotide in votes:
            allele_counts[nucleotide] += vote_weight

    # Handle "no votes at all"
    if not allele_counts:
        return None

    # Find maximum vote count
    max_count = max(allele_counts.values())

    # Get all nucleotides with maximum count and sort for deterministic output
    tied_nucleotides = sorted(
        [
            nuc
            for nuc, count in allele_counts.items()
            if math.isclose(count, max_count, rel_tol=1e-9)
        ]
    )

    if len(tied_nucleotides) > 1:
        logger.debug(
            f"Tie detected in majority voting: {tied_nucleotides} with count {max_count:.2f}"
        )

    # Return all tied nucleotides as comma-separated string
    return ",".join(tied_nucleotides)


def find_max_positive_change(frequency_changes: dict) -> Optional[str]:
    """
    Find the nucleotide(s) with the greatest positive frequency change.

    Parameters:
    -----------
    frequency_changes : dict
        Dictionary mapping nucleotide to frequency change value

    Returns:
    --------
    str or None : Nucleotide(s) with greatest positive change (comma-separated if tied),
                  or None if no positive changes exist
    """
    # Filter to only positive changes (exclude NaN and non-positive values)
    positive_changes = {
        nuc: change
        for nuc, change in frequency_changes.items()
        if not math.isnan(change) and change > 0
    }

    if not positive_changes:
        return None

    # Find maximum positive change
    max_change = max(positive_changes.values())

    # Get all nucleotides with maximum change (handle ties)
    tied_nucleotides = sorted(
        nuc
        for nuc, change in positive_changes.items()
        if math.isclose(change, max_change, rel_tol=1e-9)
    )

    if len(tied_nucleotides) > 1:
        logger.debug(
            f"Tie detected in frequency change: {tied_nucleotides} with change {max_change:.4f}"
        )

    return ",".join(tied_nucleotides)


def calculate_per_subject_frequency_changes(
    contig,
    position,
    subject_mapping,
    terminal_profile_cache,
    initial_profile_cache,
):
    """
    Calculate frequency changes per subject, then average across subjects.

    For each subjectID, computes (terminal_freq - initial_freq) for each nucleotide,
    then averages these differences across all subjectIDs with valid data.

    Parameters:
    -----------
    contig : str
        Contig identifier
    position : int
        Position on the contig
    subject_mapping : dict
        Mapping of subjectID to {'terminal': sample_id, 'initial': sample_id}
    terminal_profile_cache : dict
        Cache of terminal sample profile DataFrames indexed by (contig, position)
    initial_profile_cache : dict
        Cache of initial sample profile DataFrames indexed by (contig, position)

    Returns:
    --------
    tuple : (terminal_nucleotide, mean_frequency_changes, per_subject_changes, n_subjects_used)
        - terminal_nucleotide: str or None
        - mean_frequency_changes: dict mapping nucleotide to mean change
        - per_subject_changes: dict mapping subject to {nuc: change}
        - n_subjects_used: int
    """
    per_subject_changes = {}

    for subject, samples in subject_mapping.items():
        terminal_sample = samples["terminal"]
        initial_sample = samples["initial"]

        # Check if both samples are in cache
        if terminal_sample not in terminal_profile_cache:
            continue
        if initial_sample not in initial_profile_cache:
            continue

        terminal_profile = terminal_profile_cache[terminal_sample]
        initial_profile = initial_profile_cache[initial_sample]

        # Try to get site data from both profiles
        try:
            terminal_row = terminal_profile.loc[(contig, position)]
            initial_row = initial_profile.loc[(contig, position)]
        except KeyError:
            # Site not found in one or both profiles
            continue

        terminal_coverage = terminal_row["total_coverage"]
        initial_coverage = initial_row["total_coverage"]

        if terminal_coverage <= 0 or initial_coverage <= 0:
            continue

        # Calculate frequency change for each nucleotide for this subject
        subject_changes = {}
        for nuc in NUCLEOTIDES:
            terminal_freq = terminal_row[nuc] / terminal_coverage
            initial_freq = initial_row[nuc] / initial_coverage
            subject_changes[nuc] = terminal_freq - initial_freq

        per_subject_changes[subject] = subject_changes

    # Average changes across subjects
    n_subjects = len(per_subject_changes)
    if n_subjects == 0:
        return None, {nuc: math.nan for nuc in NUCLEOTIDES}, {}, 0

    mean_frequency_changes = {}
    for nuc in NUCLEOTIDES:
        changes = [subj_changes[nuc] for subj_changes in per_subject_changes.values()]
        mean_frequency_changes[nuc] = sum(changes) / n_subjects

    # Find nucleotide with greatest positive mean change
    terminal_nuc = find_max_positive_change(mean_frequency_changes)

    return terminal_nuc, mean_frequency_changes, per_subject_changes, n_subjects


def worker_wrapper(
    mag_data,
    profile_files,
    sample_ids,
    p_value_column,
    output_dir,
    initial_profile_files=None,
    initial_sample_ids=None,
    subject_mapping=None,
):
    """
    Worker function wrapper to unpack tuple arguments for parallel processing.

    This must be a top-level function to be picklable on macOS/Windows (spawn start method).

    Parameters:
    -----------
    mag_data : tuple
        Tuple of (mag_id, mag_sites_dataframe)
    profile_files : dict
        Dictionary mapping sample_id to profile file paths (terminal timepoint)
    sample_ids : list
        List of sample IDs to process (terminal timepoint)
    p_value_column : str
        Column name for p-values
    output_dir : str
        Output directory
    initial_profile_files : dict, optional
        Dictionary mapping sample_id to profile file paths (initial timepoint)
    initial_sample_ids : list, optional
        List of sample IDs for initial timepoint
    subject_mapping : dict, optional
        Mapping of subject_id to {'terminal': sample_id, 'initial': sample_id}

    Returns:
    --------
    dict : Results summary from process_single_mag
    """
    mag_id, mag_sites = mag_data
    return process_single_mag(
        mag_id=mag_id,
        mag_sites=mag_sites,
        profile_files=profile_files,
        sample_ids=sample_ids,
        p_value_column=p_value_column,
        output_dir=output_dir,
        initial_profile_files=initial_profile_files,
        initial_sample_ids=initial_sample_ids,
        subject_mapping=subject_mapping,
    )


# =============================================================================
# Helper functions for process_single_mag
# =============================================================================


def load_profile_cache(
    profile_files_dict: dict,
    sample_id_list: list,
    mag_id: str,
    timepoint_label: str,
) -> dict:
    """
    Load and cache profile data for a set of samples.

    Parameters:
    -----------
    profile_files_dict : dict
        Dictionary mapping sample_id to list of profile file paths
    sample_id_list : list
        List of sample IDs to load
    mag_id : str
        MAG identifier to filter profile files
    timepoint_label : str
        Label for logging (e.g., "terminal", "initial")

    Returns:
    --------
    dict : Dictionary mapping sample_id to indexed DataFrame
    """
    cache = {}
    for sample_id in sample_id_list:
        if sample_id not in profile_files_dict:
            logger.warning(
                f"No profile files found for {timepoint_label} sample {sample_id}"
            )
            continue

        # Filter profile files for this specific MAG ID
        mag_specific_files = [
            f for f in profile_files_dict[sample_id] if f"_{mag_id}_profiled.tsv" in f
        ]

        if len(mag_specific_files) == 0:
            logger.warning(
                f"No profile file found for {timepoint_label} sample {sample_id} and MAG {mag_id}"
            )
            continue
        elif len(mag_specific_files) > 1:
            raise ValueError(
                f"Multiple profile files found for {timepoint_label} sample {sample_id} "
                f"and MAG {mag_id}: {mag_specific_files}"
            )

        profile_file = mag_specific_files[0]
        profile_data = pd.read_csv(
            profile_file,
            sep="\t",
            dtype=PROFILE_DTYPES,
            usecols=PROFILE_DTYPES.keys(),
        )
        # Filter for positions with total_coverage > 0
        profile_data = profile_data[profile_data["total_coverage"] > 0]
        # Index by (contig, position) for fast lookup
        profile_data = profile_data.set_index(["contig", "position"])
        cache[sample_id] = profile_data

    return cache


def get_site_frequencies(
    contig: str,
    position: int,
    sample_ids: list,
    profile_cache: dict,
) -> tuple[dict, list, list]:
    """
    Extract nucleotide frequencies for a site across all samples.

    Parameters:
    -----------
    contig : str
        Contig identifier
    position : int
        Position on the contig
    sample_ids : list
        List of sample IDs to process
    profile_cache : dict
        Cache of sample profile DataFrames

    Returns:
    --------
    tuple : (site_frequencies, major_alleles, samples_with_site)
        - site_frequencies: dict mapping "{nuc}_{sample_id}" to frequency
        - major_alleles: list of major allele strings (or None) per sample
        - samples_with_site: list of sample_ids that have this site
    """
    site_frequencies = {}
    major_alleles = []
    samples_with_site = []

    for sample_id in sample_ids:
        if sample_id not in profile_cache:
            major_alleles.append(None)
            continue

        profile_data = profile_cache[sample_id]
        try:
            site_row = profile_data.loc[(contig, position)]
            total_coverage = site_row["total_coverage"]

            # Calculate frequencies
            sample_freqs = {nuc: site_row[nuc] / total_coverage for nuc in NUCLEOTIDES}

            # Store frequencies
            for nuc, freq in sample_freqs.items():
                site_frequencies[f"{nuc}_{sample_id}"] = freq

            # Identify major allele(s)
            max_freq = max(sample_freqs.values())
            tied_alleles = sorted(
                nuc
                for nuc, freq in sample_freqs.items()
                if math.isclose(freq, max_freq, rel_tol=1e-9)
            )
            major_alleles.append(",".join(tied_alleles))
            samples_with_site.append(sample_id)

        except KeyError:
            logger.debug(
                f"Site {contig}:{position} not found in profile for sample {sample_id}"
            )
            major_alleles.append(None)

    return site_frequencies, major_alleles, samples_with_site


def calculate_mean_frequencies(
    site_frequencies: dict,
    samples_with_site: list,
) -> dict:
    """
    Calculate mean nucleotide frequencies across samples.

    Parameters:
    -----------
    site_frequencies : dict
        Dictionary mapping "{nuc}_{sample_id}" to frequency
    samples_with_site : list
        List of sample_ids that have this site

    Returns:
    --------
    dict : Dictionary mapping nucleotide to mean frequency
    """
    if not samples_with_site:
        return {nuc: math.nan for nuc in NUCLEOTIDES}

    return {
        nuc: sum(
            site_frequencies.get(f"{nuc}_{sample_id}", 0)
            for sample_id in samples_with_site
        )
        / len(samples_with_site)
        for nuc in NUCLEOTIDES
    }


def determine_terminal_by_mean_freq(mean_frequencies: dict) -> Optional[str]:
    """
    Determine terminal nucleotide based on highest mean frequency.

    Parameters:
    -----------
    mean_frequencies : dict
        Dictionary mapping nucleotide to mean frequency

    Returns:
    --------
    str or None : Terminal nucleotide(s), comma-separated if tied
    """
    # Check for all NaN
    valid_freqs = {
        nuc: freq for nuc, freq in mean_frequencies.items() if not math.isnan(freq)
    }
    if not valid_freqs:
        return None

    max_freq = max(valid_freqs.values())
    tied = sorted(
        nuc
        for nuc, freq in valid_freqs.items()
        if math.isclose(freq, max_freq, rel_tol=1e-9)
    )
    return ",".join(tied) if tied else None


def split_count_nucleotides(series: pd.Series) -> dict:
    """
    Count nucleotides in a series, splitting comma-separated ties.

    Parameters:
    -----------
    series : pd.Series
        Series of nucleotide strings (possibly comma-separated)

    Returns:
    --------
    dict : Dictionary mapping nucleotide to count
    """
    counter: dict[str, int] = {}
    for val in series.dropna():
        for nuc in str(val).split(","):
            if nuc and nuc in NUCLEOTIDES:
                counter[nuc] = counter.get(nuc, 0) + 1
    return counter


def process_single_mag(
    mag_id,
    mag_sites,
    profile_files,
    sample_ids,
    p_value_column,
    output_dir,
    initial_profile_files=None,
    initial_sample_ids=None,
    subject_mapping=None,
):
    """
    Process a single MAG's significant sites.

    Parameters:
    -----------
    mag_id : str
        MAG identifier
    mag_sites : pd.DataFrame
        Sites for this MAG
    profile_files : dict
        Dictionary mapping sample_id to profile file paths (terminal timepoint)
    sample_ids : list
        List of sample IDs to process (terminal timepoint)
    p_value_column : str
        Column name for p-values
    output_dir : str
        Output directory
    initial_profile_files : dict, optional
        Dictionary mapping sample_id to profile file paths (initial timepoint)
    initial_sample_ids : list, optional
        List of sample IDs for initial timepoint
    subject_mapping : dict, optional
        Mapping of subject_id to {'terminal': sample_id, 'initial': sample_id}

    Returns:
    --------
    dict : Results summary
    """
    use_freq_change = subject_mapping is not None and len(subject_mapping) > 0

    # Load profile caches
    terminal_cache = load_profile_cache(profile_files, sample_ids, mag_id, "terminal")

    initial_cache = {}
    if use_freq_change and initial_profile_files and initial_sample_ids:
        initial_cache = load_profile_cache(
            initial_profile_files, initial_sample_ids, mag_id, "initial"
        )

    # Process each significant site
    results_mean_freq = []
    results_majority_vote = []
    results_freq_change = []
    intermediate_data = []

    for _, site in mag_sites.iterrows():
        contig = site["contig"]
        position = site["position"]

        # Get terminal frequencies
        site_freqs, major_alleles, samples_with_site = get_site_frequencies(
            contig, position, sample_ids, terminal_cache
        )

        mean_freqs = calculate_mean_frequencies(site_freqs, samples_with_site)

        # Method 1: Mean frequency
        terminal_mean_freq = determine_terminal_by_mean_freq(mean_freqs)

        # Method 2: Majority voting
        terminal_majority = perform_majority_voting(major_alleles)

        # Method 3: Frequency change (per-subject, then averaged)
        terminal_freq_change = None
        freq_changes = {nuc: math.nan for nuc in NUCLEOTIDES}
        n_subjects_used = 0

        if use_freq_change:
            (
                terminal_freq_change,
                freq_changes,
                _per_subj,
                n_subjects_used,
            ) = calculate_per_subject_frequency_changes(
                contig, position, subject_mapping, terminal_cache, initial_cache
            )

        # Build intermediate data row
        site_row_data = {
            "mag_id": mag_id,
            "contig": contig,
            "position": position,
            "gene_id": site["gene_id"],
            p_value_column: site[p_value_column],
            # Terminal sample frequencies
            **{
                f"{nuc}_{sid}": site_freqs.get(f"{nuc}_{sid}", math.nan)
                for nuc in NUCLEOTIDES
                for sid in sample_ids
            },
            # Mean frequencies (terminal)
            **{f"{nuc}_mean_frequency": mean_freqs[nuc] for nuc in NUCLEOTIDES},
            "n_samples_used_for_mean": len(samples_with_site),
            # Major alleles per sample
            **{
                f"major_allele_{sample_ids[i]}": major_alleles[i]
                for i in range(len(sample_ids))
                if i < len(major_alleles)
            },
            "terminal_nucleotide_mean_freq": terminal_mean_freq,
            "terminal_nucleotide_majority_vote": terminal_majority,
        }

        if use_freq_change:
            site_row_data.update(
                {
                    **{
                        f"{nuc}_frequency_change": freq_changes[nuc]
                        for nuc in NUCLEOTIDES
                    },
                    "n_subjects_used_for_freq_change": n_subjects_used,
                    "terminal_nucleotide_freq_change": terminal_freq_change,
                }
            )

        intermediate_data.append(site_row_data)
        results_mean_freq.append(terminal_mean_freq)
        results_majority_vote.append(terminal_majority)
        results_freq_change.append(terminal_freq_change)

    # Build output DataFrames
    processed_sites = mag_sites.copy()
    processed_sites["terminal_nucleotide_mean_freq"] = results_mean_freq
    processed_sites["terminal_nucleotide_majority_vote"] = results_majority_vote

    if use_freq_change:
        processed_sites["terminal_nucleotide_freq_change"] = results_freq_change

    if "source_file" in processed_sites.columns:
        processed_sites = processed_sites.drop(columns=["source_file"])

    intermediate_df = pd.DataFrame(intermediate_data)

    # Save outputs
    mag_output_dir = Path(output_dir) / mag_id
    mag_output_dir.mkdir(parents=True, exist_ok=True)

    main_output_file = mag_output_dir / f"{mag_id}_terminal_nucleotides.tsv"
    processed_sites.to_csv(main_output_file, sep="\t", index=False)

    intermediate_output_file = (
        mag_output_dir / f"{mag_id}_nucleotide_frequencies.tsv.gz"
    )
    intermediate_df.to_csv(
        intermediate_output_file, sep="\t", index=False, compression="gzip"
    )

    # Summary statistics
    terminal_counts = split_count_nucleotides(
        processed_sites["terminal_nucleotide_mean_freq"]
    )
    majority_counts = split_count_nucleotides(
        processed_sites["terminal_nucleotide_majority_vote"]
    )

    freq_change_counts = {}
    if use_freq_change:
        freq_change_counts = split_count_nucleotides(
            processed_sites["terminal_nucleotide_freq_change"]
        )

    result = {
        "mag_id": mag_id,
        "sites_processed": len(processed_sites),
        "mean_freq_terminal_nucleotides": terminal_counts,
        "majority_vote_terminal_nucleotides": majority_counts,
        "main_output": str(main_output_file),
        "intermediate_output": (
            str(intermediate_output_file) if not intermediate_df.empty else None
        ),
        "samples_processed": len(sample_ids),
    }

    if use_freq_change:
        result["freq_change_terminal_nucleotides"] = freq_change_counts
        result["subjects_used"] = len(subject_mapping)

    return result


def main():
    """Main function for terminal nucleotide analysis CLI."""
    # Parse arguments first to get log level
    parser = argparse.ArgumentParser(
        description="Identify terminal nucleotides in specified group and timepoint samples using dual methodology",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--significant_sites",
        help="Tab-separated file containing significant sites with columns 'mag_id', 'contig', 'position', 'gene_id', and either 'min_p_value' or 'q_value'",
        type=str,
        required=True,
        metavar="filepath",
    )

    parser.add_argument(
        "--profile_dir",
        help="Directory containing nucleotide frequency profiles for each sample (organized in subdirectories by sample_id)",
        type=str,
        required=True,
        metavar="directory",
    )

    parser.add_argument(
        "--group",
        help="Specific group name to identify terminal timepoint samples",
        type=str,
        required=True,
        metavar="string",
    )

    parser.add_argument(
        "--timepoint",
        help="Specific timepoint to identify terminal samples",
        type=str,
        required=True,
        metavar="string",
    )

    parser.add_argument(
        "--initial_timepoint",
        help="Optional initial timepoint for frequency change analysis. When provided, enables "
        "calculation of terminal nucleotide based on greatest positive frequency change "
        "(terminal - initial). Samples from the same group at this timepoint will be used.",
        type=str,
        default=None,
        metavar="string",
    )

    parser.add_argument(
        "--metadata",
        help="Path to metadata table containing sample information",
        type=str,
        required=True,
        metavar="filepath",
    )

    parser.add_argument(
        "--output",
        help="Output directory for results",
        type=str,
        required=True,
        metavar="directory",
    )

    parser.add_argument(
        "--p_value_column",
        help="Column to use for p-value filtering",
        type=str,
        choices=["min_p_value", "q_value"],
        default="q_value",
        metavar="string",
    )

    parser.add_argument(
        "--p_value_threshold",
        type=float,
        default=0.05,
        help="Threshold for significance filtering using the selected p-value column",
        metavar="float",
    )

    parser.add_argument(
        "--test-type",
        type=str,
        choices=[
            "two_sample_unpaired_tTest",
            "two_sample_unpaired_MannWhitney",
            "two_sample_unpaired_tTest_abs",
            "two_sample_unpaired_MannWhitney_abs",
            "two_sample_paired_tTest",
            "two_sample_paired_Wilcoxon",
            "two_sample_paired_tTest_abs",
            "two_sample_paired_Wilcoxon_abs",
            "single_sample_tTest",
            "single_sample_Wilcoxon",
            "CMH",
            "LMM",
            "LMM_abs",
            # "lmm_across_time",
            # "cmh_across_time",
        ],
        default="two_sample_paired_tTest",
        help="A test type to filter the `significant_sites` dataframe.",
    )

    parser.add_argument(
        "--group_analyzed",
        help="Optional group name to filter significant sites by the 'group_analyzed' column if present",
        type=str,
        default=None,
        metavar="string",
    )

    parser.add_argument(
        "--cpus",
        help="Number of CPU cores to use for parallel processing (default: all available)",
        type=int,
        default=multiprocessing.cpu_count(),
        metavar="integer",
    )

    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level",
        metavar="string",
    )

    args = parser.parse_args()

    # Setup logging with the specified level
    setup_logging(level=getattr(logging, args.log_level))

    # Validate input files exist
    input_files = [args.significant_sites, args.metadata]
    for file_path in input_files:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")

    if not Path(args.profile_dir).exists():
        raise FileNotFoundError(f"Profile directory not found: {args.profile_dir}")

    # Create output directory
    output_path = Path(args.output)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info("Starting enhanced terminal nucleotide analysis")

    # Load significant sites
    logger.info("Loading significant sites table")
    significant_sites = pd.read_csv(args.significant_sites, sep="\t")

    if significant_sites.empty:
        raise ValueError("Significant sites file is empty")

    # Validate required columns
    required_cols = ["mag_id", "contig", "position", "gene_id", args.p_value_column]
    missing_cols = set(required_cols) - set(significant_sites.columns)
    if missing_cols:
        raise ValueError(
            f"Missing required columns in significant sites: {missing_cols}"
        )

    # Apply significance threshold filtering
    logger.info(
        f"Applying significance threshold: {args.p_value_column} <= {args.p_value_threshold}"
    )
    initial_count = len(significant_sites)
    significant_sites = significant_sites[
        significant_sites[args.p_value_column] <= args.p_value_threshold
    ].copy()
    logger.info(
        f"After significance filtering ({args.p_value_column} <= {args.p_value_threshold}), "
        f"{len(significant_sites):,}/{initial_count:,} sites remain"
    )

    if significant_sites.empty:
        raise ValueError(
            f"No sites remain after applying {args.p_value_column} <= {args.p_value_threshold} threshold"
        )

    # Filter by test_type if specified
    if args.test_type:
        if "test_type" in significant_sites.columns:
            initial_count = len(significant_sites)
            significant_sites = significant_sites[
                significant_sites["test_type"] == args.test_type
            ].copy()
            logger.info(
                f"Filtered by test_type '{args.test_type}': {len(significant_sites):,}/{initial_count:,} sites remain"
            )
            if significant_sites.empty:
                raise ValueError(f"No sites found for test_type '{args.test_type}'")
        else:
            logger.warning(
                f"test_type column not found in significant_sites file. Ignoring --test_type filter."
            )

    # Filter by group_analyzed if specified
    if args.group_analyzed:
        if "group_analyzed" in significant_sites.columns:
            initial_count = len(significant_sites)
            significant_sites = significant_sites[
                significant_sites["group_analyzed"] == args.group_analyzed
            ].copy()
            logger.info(
                f"Filtered by group_analyzed '{args.group_analyzed}': {len(significant_sites)}/{initial_count} sites remain"
            )
            if significant_sites.empty:
                raise ValueError(
                    f"No sites found for group_analyzed '{args.group_analyzed}'"
                )
        else:
            logger.warning(
                f"group_analyzed column not found in significant_sites file. Ignoring --group_analyzed filter."
            )

    # Get matching samples from metadata
    terminal_samples = get_terminal_samples(args.metadata, args.group, args.timepoint)

    # Get initial timepoint samples and subject mapping if frequency change analysis is requested
    initial_samples = None
    initial_profile_files = None
    subject_mapping = None
    use_freq_change = args.initial_timepoint is not None

    if use_freq_change:
        logger.info(
            f"Frequency change analysis enabled: comparing timepoint '{args.timepoint}' vs '{args.initial_timepoint}'"
        )
        # Get subject-to-sample mapping for paired analysis
        subject_mapping = get_subject_sample_mapping(
            args.metadata, args.group, args.timepoint, args.initial_timepoint
        )

        # Extract initial sample IDs from subject mapping
        initial_samples = [mapping["initial"] for mapping in subject_mapping.values()]
        logger.info(
            f"Found {len(subject_mapping)} subject(s) with paired samples for frequency change analysis"
        )

    # Find profile files for specific MAG IDs from significant sites
    mag_groups = significant_sites.groupby("mag_id")
    mag_ids = list(mag_groups.groups.keys())
    profile_files = find_profile_files(args.profile_dir, terminal_samples, mag_ids)

    if not profile_files:
        raise ValueError(
            f"No profile files found for terminal samples in {args.profile_dir}"
        )

    # Find profile files for initial timepoint samples
    if use_freq_change:
        initial_profile_files = find_profile_files(
            args.profile_dir, initial_samples, mag_ids
        )
        if not initial_profile_files:
            raise ValueError(
                f"No profile files found for initial timepoint samples in {args.profile_dir}"
            )

    # Optimize CPU usage: use minimum of specified CPUs and number of MAGs
    # to avoid wasting resources when there are fewer MAGs than CPUs
    effective_cpus = min(args.cpus, len(mag_ids))

    logger.info(
        f"Processing {len(mag_ids)} MAG(s) with significant sites using {effective_cpus} CPU(s)"
    )

    # Prepare MAG data for parallel processing
    # Create list of tuples: [(mag_id, mag_sites_dataframe), (mag_id, mag_sites_dataframe), ...]
    # Each tuple contains the MAG ID and its corresponding DataFrame of significant sites
    mag_data_list = [(mag_id, mag_groups.get_group(mag_id)) for mag_id in mag_ids]

    # Process MAGs in parallel using multiprocessing
    all_results = []
    with multiprocessing.Pool(processes=effective_cpus) as pool:
        # Use functools.partial to bind shared arguments to top-level worker function
        # This ensures the worker is picklable on macOS/Windows (spawn start method)
        func = partial(
            worker_wrapper,
            profile_files=profile_files,
            sample_ids=terminal_samples,
            p_value_column=args.p_value_column,
            output_dir=args.output,
            initial_profile_files=initial_profile_files,
            initial_sample_ids=initial_samples,
            subject_mapping=subject_mapping,
        )

        # Execute parallel processing with streaming results using imap_unordered
        # Results are yielded as they complete rather than waiting for all to finish
        for result in tqdm(
            pool.imap_unordered(func, mag_data_list),
            total=len(mag_data_list),
            desc="Processing MAGs",
        ):
            # Filter out None results (failed MAGs) and collect successful ones
            if result is not None:
                all_results.append(result)

    # Create summary report
    if all_results:
        summary_data = []
        for result in all_results:
            row_data = {
                "mag_id": result["mag_id"],
                "sites_processed": result["sites_processed"],
                "samples_processed": result["samples_processed"],
                "main_output": result["main_output"],
                "intermediate_output": result["intermediate_output"],
                **{
                    f"mean_freq_terminal_{nuc}_count": count
                    for nuc, count in result["mean_freq_terminal_nucleotides"].items()
                },
                **{
                    f"majority_vote_terminal_{nuc}_count": count
                    for nuc, count in result[
                        "majority_vote_terminal_nucleotides"
                    ].items()
                },
            }

            # Add frequency change results if available
            if "freq_change_terminal_nucleotides" in result:
                row_data["subjects_used"] = result.get("subjects_used", 0)
                row_data.update(
                    {
                        f"freq_change_terminal_{nuc}_count": count
                        for nuc, count in result[
                            "freq_change_terminal_nucleotides"
                        ].items()
                    }
                )

            summary_data.append(row_data)

        summary_df = pd.DataFrame(summary_data)
        summary_file = output_path / "terminal_nucleotide_analysis_summary.tsv"
        summary_df.to_csv(summary_file, sep="\t", index=False)

        logger.info(f"Analysis complete. Processed {len(all_results)} MAG(s)")
        logger.info(f"Summary report saved to: {summary_file}")

        # Log overall statistics
        total_sites = sum(r["sites_processed"] for r in all_results)
        total_mags = len(all_results)
        logger.info(f"Total: {total_sites:,} site(s) across {total_mags} MAG(s)")

        if use_freq_change:
            logger.info(
                "All three methods completed: mean frequency, majority voting, and frequency change"
            )
        else:
            logger.info(
                "Both mean frequency and majority voting methods completed successfully"
            )

    else:
        logger.warning("No MAGs were successfully processed")


if __name__ == "__main__":
    main()
