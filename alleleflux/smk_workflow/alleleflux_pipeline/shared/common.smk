"""Common Snakemake utilities and configuration.

This module provides shared configuration, helper functions, and wildcard constraints
for the AlleleFlux Snakemake pipeline. It handles:
- Global configuration parsing (data type, output directory, timepoints, groups)
- Resource management (memory parsing, per-rule overrides)
- MAG eligibility functions for QC and preprocessing stages
- Sample metadata parsing for longitudinal analysis
- Wildcard constraints for consistent rule matching
"""

import os
import pandas as pd
from glob import glob
from collections import defaultdict
from pathlib import Path
from snakemake.logging import logger

# Load the configuration file
# configfile: os.path.join(workflow.basedir, "config.yml")


# =============================================================================
# Resource Management
# =============================================================================

def parse_mem(mem_value):
    """
    Convert memory string to MB for Snakemake resources.
    
    Supports formats: "8G", "8GB", "8192M", "8192MB", "8192" (assumes MB)
    Case-insensitive. Uses binary units (1 GB = 1024 MB).
    
    Args:
        mem_value: Memory value as string (e.g., "8G") or int (MB)
    
    Returns:
        int: Memory in MB
    
    Examples:
        >>> parse_mem("8G")
        8192
        >>> parse_mem("100GB")
        102400
        >>> parse_mem("8192M")
        8192
        >>> parse_mem(8192)
        8192
    """
    if isinstance(mem_value, (int, float)):
        return int(mem_value)
    
    mem_str = str(mem_value).strip().upper()
    
    # Remove 'B' suffix if present (e.g., "8GB" -> "8G")
    if mem_str.endswith("B"):
        mem_str = mem_str[:-1]
    
    if mem_str.endswith("G"):
        return int(float(mem_str[:-1]) * 1024)
    elif mem_str.endswith("M"):
        return int(float(mem_str[:-1]))
    elif mem_str.endswith("K"):
        return max(1, int(float(mem_str[:-1]) / 1024))
    else:
        # Assume MB if no unit
        return int(float(mem_str))


def parse_time_to_minutes(time_str):
    """
    Convert a time string (HH:MM:SS or D-HH:MM:SS) to total minutes.
    
    Args:
        time_str: Time string in HH:MM:SS or D-HH:MM:SS format
    
    Returns:
        int: Total minutes
    
    Examples:
        >>> parse_time_to_minutes("24:00:00")
        1440
        >>> parse_time_to_minutes("4:00:00")
        240
        >>> parse_time_to_minutes("1-00:00:00")
        1440
        >>> parse_time_to_minutes("0:30:00")
        30
    """
    time_str = str(time_str).strip()
    
    # Handle D-HH:MM:SS format
    days = 0
    if "-" in time_str:
        day_part, time_str = time_str.split("-", 1)
        days = int(day_part)
    
    parts = time_str.split(":")
    if len(parts) == 3:
        hours, minutes, seconds = int(parts[0]), int(parts[1]), int(parts[2])
    elif len(parts) == 2:
        hours, minutes = 0, int(parts[0])
        seconds = int(parts[1])
    else:
        raise ValueError(f"Invalid time format: {time_str}. Use HH:MM:SS or D-HH:MM:SS.")
    
    return days * 24 * 60 + hours * 60 + minutes + (1 if seconds > 0 else 0)


def minutes_to_time_str(total_minutes):
    """
    Convert total minutes to HH:MM:SS or D-HH:MM:SS format.
    
    Uses D-HH:MM:SS format only when days >= 1 for cleaner output.
    
    Args:
        total_minutes: Total minutes as int
    
    Returns:
        str: Time string in HH:MM:SS or D-HH:MM:SS format
    
    Examples:
        >>> minutes_to_time_str(240)
        '4:00:00'
        >>> minutes_to_time_str(1440)
        '1-00:00:00'
        >>> minutes_to_time_str(1500)
        '1-01:00:00'
    """
    total_minutes = max(1, int(total_minutes))
    days = total_minutes // (24 * 60)
    remaining = total_minutes % (24 * 60)
    hours = remaining // 60
    minutes = remaining % 60
    
    if days > 0:
        return f"{days}-{hours:02d}:{minutes:02d}:00"
    else:
        return f"{hours}:{minutes:02d}:00"


def get_resource(rule_name, resource_type, default=None):
    """
    Get resource value for a rule, checking for per-rule override first.
    
    This enables the 'escape hatch' pattern where power users can override
    resources for specific rules via the resources_override config section,
    while most users just use the flat defaults.
    
    Args:
        rule_name: Name of the Snakemake rule (e.g., "profile", "qc", "statistical_tests")
        resource_type: Type of resource ("threads_per_job", "mem_per_job", "time")
        default: Default value if not specified anywhere (uses config default if None)
    
    Returns:
        Resource value (parsed if memory)
    
    Example:
        # In a rule:
        resources:
            mem_mb=get_resource("profile", "mem_per_job"),
            time=get_resource("profile", "time")
    """
    # Check for per-rule override first
    overrides = config.get("resources_override", {}).get(rule_name, {})
    if resource_type in overrides:
        value = overrides[resource_type]
    elif default is not None:
        value = default
    else:
        # Get from main resources section
        value = config.get("resources", {}).get(resource_type)
    
    # Parse memory values to MB
    if resource_type == "mem_per_job" and value is not None:
        return parse_mem(value)
    
    return value


# -- Internal helpers for retry / resource stepping --

def _get_retries_for_rule(rule_name=None):
    """Get retry count, checking per-rule override first, then global default."""
    if rule_name:
        overrides = config.get("resources_override", {}).get(rule_name, {})
        if "retries" in overrides:
            return int(overrides["retries"])
    return int(config.get("resources", {}).get("retries", 2))


def _get_step_value(rule_name, key, global_default):
    """Get a step value (mem_step, time_step) with per-rule override support."""
    if rule_name:
        overrides = config.get("resources_override", {}).get(rule_name, {})
        if key in overrides:
            return overrides[key]
    return config.get("resources", {}).get(key, global_default)


# -- Public resource functions used by rules --

def get_threads(rule_name=None):
    """Get threads_per_job, optionally checking for rule-specific override."""
    if rule_name:
        return get_resource(rule_name, "threads_per_job")
    return config.get("resources", {}).get("threads_per_job", 1)


def get_retries(rule_name=None):
    """
    Get retry count for a rule, for use in the rule-level retries: directive.
    
    Checks per-rule override first, then falls back to global resources.retries.
    Default is 2 retries if not specified.
    
    Args:
        rule_name: Name of the Snakemake rule (optional)
    
    Returns:
        int: Number of retries
    """
    return _get_retries_for_rule(rule_name)


def get_mem_mb(rule_name=None):
    """
    Get memory in MB, with automatic scaling on retry when mem_step is configured.
    
    When retries > 0 and mem_step is set, returns a callable
    ``lambda wildcards, attempt: base + (attempt-1) * step`` so that
    Snakemake allocates more memory on each retry.
    When no stepping is configured, returns a static int (backward-compatible).
    
    Args:
        rule_name: Name of the Snakemake rule (optional)
    
    Returns:
        int or callable: Memory in MB (static or attempt-aware)
    """
    if rule_name:
        base_mb = get_resource(rule_name, "mem_per_job")
    else:
        base_mb = parse_mem(config.get("resources", {}).get("mem_per_job", "8G"))
    
    retries = _get_retries_for_rule(rule_name)
    step_raw = _get_step_value(rule_name, "mem_step", None)
    
    if retries > 0 and step_raw:
        step_mb = parse_mem(step_raw)
        if step_mb > 0:
            return lambda wildcards, attempt: base_mb + (attempt - 1) * step_mb
    
    return base_mb


def get_time(rule_name=None):
    """
    Get wall time, with automatic scaling on retry when time_step is configured.
    
    When retries > 0 and time_step is set, returns a callable
    ``lambda wildcards, attempt: base_time + (attempt-1) * step`` so that
    Snakemake allocates more time on each retry.
    When no stepping is configured, returns a static string (backward-compatible).
    
    Args:
        rule_name: Name of the Snakemake rule (optional)
    
    Returns:
        str or callable: Wall time (static or attempt-aware)
    """
    if rule_name:
        base_time = get_resource(rule_name, "time")
    else:
        base_time = config.get("resources", {}).get("time", "24:00:00")
    
    retries = _get_retries_for_rule(rule_name)
    step_raw = _get_step_value(rule_name, "time_step", None)
    
    if retries > 0 and step_raw:
        base_mins = parse_time_to_minutes(base_time)
        step_mins = parse_time_to_minutes(step_raw)
        if step_mins > 0:
            return lambda wildcards, attempt: minutes_to_time_str(
                base_mins + (attempt - 1) * step_mins
            )
    
    return base_time


# =============================================================================
# Global Constants
# =============================================================================

# Taxonomy levels for aggregation (order matters - from broad to specific)
TAXONOMY_LEVELS = ["phylum", "class", "order", "family", "genus", "species"]

# Statistical test types - centralized for consistency
BETWEEN_GROUP_TEST_TYPES = ["two_sample_unpaired", "two_sample_paired", "lmm", "cmh"]
WITHIN_GROUP_TEST_TYPES = ["single_sample", "lmm_across_time", "cmh_across_time"]
ALL_TEST_TYPES = BETWEEN_GROUP_TEST_TYPES + WITHIN_GROUP_TEST_TYPES

# =============================================================================
# Global Configuration
# =============================================================================

# Define the global data_type variable to be used across all Snakemake files
DATA_TYPE = config["analysis"]["data_type"]
OUTDIR = config["output"]["root_dir"]
DN_DS_TEST_TYPE = config["dnds"]["dn_ds_test_type"]


def get_base_test_type(test_type):
    """
    Convert a specific test type to its base eligibility test type.
    
    Maps specific test variants (e.g., 'two_sample_paired_tTest') to their
    base types (e.g., 'two_sample_paired') for eligibility checking.
    
    Parameters:
        test_type: The specific test type string from config
    
    Returns:
        The base test type string for eligibility lookup
    
    Raises:
        ValueError: If the test type is not recognized
    """
    if test_type in ["two_sample_unpaired_tTest", "two_sample_unpaired_MannWhitney", 
                     "two_sample_unpaired_tTest_abs", "two_sample_unpaired_MannWhitney_abs"]:
        return "two_sample_unpaired"
    elif test_type in ["two_sample_paired_tTest", "two_sample_paired_Wilcoxon", 
                       "two_sample_paired_tTest_abs", "two_sample_paired_Wilcoxon_abs"]:
        return "two_sample_paired"
    elif test_type in ["single_sample_tTest", "single_sample_Wilcoxon"]:
        return "single_sample"
    elif test_type in ["lmm", "lmm_abs", "lmm_across_time", "cmh", "cmh_across_time"]:
        return test_type
    else:
        raise ValueError(f"Unsupported test type: {test_type}")


if DATA_TYPE == "single":
    OUTDIR = os.path.join(OUTDIR, "single_timepoint")
elif DATA_TYPE == "longitudinal":
    OUTDIR = os.path.join(OUTDIR, "longitudinal")

# Profile reuse configuration
# If profiles_path is specified and exists, use existing profiles instead of generating new ones
EXISTING_PROFILES_PATH = config["input"].get("profiles_path", "")
USE_EXISTING_PROFILES = bool(EXISTING_PROFILES_PATH and os.path.isdir(EXISTING_PROFILES_PATH))
PROFILES_DIR = EXISTING_PROFILES_PATH if USE_EXISTING_PROFILES else os.path.join(OUTDIR, "profiles")

if USE_EXISTING_PROFILES:
    logger.info(f"Using existing profiles from: {EXISTING_PROFILES_PATH}")
else:
    logger.info(f"Profiles will be generated in: {PROFILES_DIR}")

timepoints_labels = []
focus_timepoints = {}

for time_combo in config["analysis"]["timepoints_combinations"]:
    # Standardized format: all entries are dictionaries with "timepoint" key
    timepoint = time_combo["timepoint"]
    
    if len(timepoint) == 1 and DATA_TYPE == "single":
        # Single timepoint
        tp = timepoint[0]
        timepoints_labels.append(tp)
        # For single data type, the timepoint itself is the focus
        focus_timepoints[tp] = tp
    elif len(timepoint) == 2 and DATA_TYPE == "longitudinal":
        # Multiple timepoints (for longitudinal analysis)
        label = f"{timepoint[0]}_{timepoint[1]}"
        timepoints_labels.append(label)
        if "focus" in time_combo:
            # Ensure focus is one of the two timepoints
            if time_combo["focus"] in timepoint:
                focus_timepoints[label] = time_combo["focus"]
            else:
                raise ValueError(f"Invalid focus timepoint '{time_combo['focus']}' for combination '{label}'. Must be one of {timepoint}")
        else:
            # Default to second timepoint if focus not specified
            logger.warning(f"No focus specified for timepoint combination {label}. Using '{timepoint[1]}' as default.")
            focus_timepoints[label] = timepoint[1]

# Define valid focus timepoint values for each timepoint label
valid_focus_timepoints = {}
for tp in timepoints_labels:
    if "_" in tp and DATA_TYPE == "longitudinal":  # It's a two-timepoint combination with focus
        valid_focus_timepoints[tp] = tp.split("_")
    elif DATA_TYPE == "single":
        # For single timepoint, the timepoint itself is the focus
        valid_focus_timepoints[tp] = [tp]

groups_labels = [
    f"{gr['treatment']}_{gr['control']}" for gr in config["analysis"]["groups_combinations"]
]  # ["G1_G2", "G2_G4"]

# Flatten the list of group values from the groups_combinations config.
group_values = sorted(
    {gr["treatment"] for gr in config["analysis"]["groups_combinations"]}
    | {gr["control"] for gr in config["analysis"]["groups_combinations"]}
)
# Build the regex: optionally an underscore and one of the allowed group values.
group_str_regex = "(_({}))?".format("|".join(group_values))

# =============================================================================
# Per-timepoint cache helpers
# =============================================================================
# The ``compute_allele_freq_per_timepoint`` rule writes one Parquet cache file
# per (MAG, gr_combo, single timepoint).  Scoping the cache to gr_combo is
# required by the checkpoint architecture: each (tp_combo, gr_combo) resolves
# its eligibility_table checkpoint separately, so only the QC files for ONE
# combination are guaranteed to exist when the cache is first needed.  By
# keying the cache on gr_combo we can use a SINGLE canonical QC file (the
# first tp_combo in config order that contains this timepoint) as the input,
# which Snakemake resolves as a normal rule dependency — no cross-combination
# QC gathering needed.
#
# Deduplication is within each gr_combo, across tp_combos: all tp_combos that
# share the same (gr_combo, timepoint) reuse the same Parquet cache.  Cache files
# are NOT shared across gr_combos — each gr_combo builds its own set.
# For drido (15 tp_combos, 6 gr_combos, 8 unique timepoints):
#   Before: 15 tp_combos × 6 gr_combos × 2 timepoints = 180 profile-read operations
#   After:  8 unique timepoints × 6 gr_combos          =  48 cache-write jobs

# unique_timepoints: ordered list of individual timepoints, each appearing exactly once.
#
# timepoints_labels contains combo labels like "5mo_10mo", "5mo_16mo", "8mo_10mo".
# This loop splits each label into its constituent timepoints ("5mo", "10mo") and
# collects them in first-appearance order using _seen_tps as a visited set.
#
# Example (drido config, 15 tp_combos):
#   "5mo_10mo" → adds "5mo", "10mo"
#   "5mo_16mo" → "5mo" already seen; adds "16mo"
#   "8mo_10mo" → "10mo" already seen; adds "8mo"   ← 8mo appears last (Strategy C)
#   Result: ["5mo", "10mo", "16mo", "22mo", "28mo", "34mo", "40mo", "8mo"]
#
# This list becomes the {timepoint} wildcard constraint — the single-timepoint
# wildcard used in cache filenames, distinct from {timepoints} (the combo label).
unique_timepoints = []
_seen_tps = set()
for tp_label in timepoints_labels:
    if DATA_TYPE == "longitudinal":
        parts = tp_label.split("_")  # "5mo_10mo" → ["5mo", "10mo"]
    else:
        parts = [tp_label]           # single data: label is already one timepoint
    for tp in parts:
        if tp not in _seen_tps:
            _seen_tps.add(tp)
            unique_timepoints.append(tp)

# timepoint_gr_to_canonical_tp: maps (individual_timepoint, gr_combo) → the FIRST
# tp_combo label in config order that contains that timepoint.
#
# WHAT THIS ACHIEVES (plain language):
#   For each group combination, we build one Parquet cache file per timepoint.
#   Profile TSVs for that timepoint are read exactly once per gr_combo — not once
#   per tp_combo that happens to include it.  Every tp_combo that shares the same
#   (gr_combo, timepoint) reads from the same cache instead of re-reading raw profiles.
#
#   Cache files are NOT shared across gr_combos.
#   MAG1_1D_AL_5mo_allele_frequency.parquet  ← built once for 1D_AL
#   MAG1_2D_AL_5mo_allele_frequency.parquet  ← built once for 2D_AL, separately
#   These are independent because each gr_combo compares different sets of samples.
#
#   The deduplication is within each gr_combo, across tp_combos:
#   drido example — "5mo" appears in 6 tp_combos (5mo_10mo, 5mo_16mo, …, 5mo_40mo).
#   Without caching: 5mo profiles read 6 tp_combos × 6 gr_combos = 36 times per MAG.
#   With caching:    5mo profiles read 1 time per gr_combo   × 6 gr_combos = 6 times.
#
# WHY THIS EXISTS — the checkpoint problem:
#   Snakemake resolves the eligibility_table checkpoint one (tp_combo, gr_combo) at a
#   time.  When the cache rule for "5mo" fires during the "5mo_10mo-1D_AL" combination,
#   only QC_5mo_10mo-1D_AL/ exists on disk.  QC_5mo_16mo-1D_AL/ may not yet exist.
#   The cache rule therefore cannot depend on QC files from other combinations.
#
#   The fix: use ONE canonical QC file per (gr_combo, timepoint) — the file from the
#   first tp_combo (in config order) that contains that timepoint.  This file is built
#   by a regular rule chain (qc → generate_metadata), not gated by any other
#   combination's checkpoint, so it is always resolvable.
#
# HOW IT IS BUILT:
#   Walk timepoints_labels in config order.  For each tp_label (e.g. "5mo_10mo"),
#   split into individual timepoints.  For each (tp, gr_combo) pair not yet in the
#   dict, record tp_label as its canonical source.  The guard
#   `if key not in timepoint_gr_to_canonical_tp` means later tp_labels that also
#   contain the same timepoint are ignored — first-seen wins.
#
# Example (drido, showing 1D_AL only — all 6 gr_combos behave identically):
#   tp_label="5mo_10mo": ("5mo","1D_AL")→"5mo_10mo", ("10mo","1D_AL")→"5mo_10mo"
#   tp_label="5mo_16mo": ("5mo","1D_AL") already set; ("16mo","1D_AL")→"5mo_16mo"
#   tp_label="10mo_16mo": both already set → nothing added
#   tp_label="8mo_10mo": ("8mo","1D_AL")→"8mo_10mo";  ("10mo","1D_AL") already set
#
# The resulting entry for ("5mo","1D_AL") is always "5mo_10mo" regardless of how
# many later tp_combos also contain "5mo".
#
# See docs/source/usage/allele_freq_cache_architecture.md for the full walkthrough.
timepoint_gr_to_canonical_tp = {}
for tp_label in timepoints_labels:
    if DATA_TYPE == "longitudinal":
        parts = tp_label.split("_")  # "5mo_10mo" → ["5mo", "10mo"]
    else:
        parts = [tp_label]
    for tp in parts:
        for gr_label in groups_labels:
            key = (tp, gr_label)
            if key not in timepoint_gr_to_canonical_tp:
                timepoint_gr_to_canonical_tp[key] = tp_label  # first in config order wins

# =============================================================================
# Wildcard Constraints
# =============================================================================

wildcard_constraints:
    groups=f"({'|'.join(groups_labels)})",
    timepoints=f"({'|'.join(timepoints_labels)})",
    timepoint=f"({'|'.join(unique_timepoints)})",
    taxon="(" + "|".join(TAXONOMY_LEVELS) + "|domain)",
    test_type="(" + "|".join(ALL_TEST_TYPES) + "|)",
    sub_test="(MannWhitney|Wilcoxon|tTest|LMM|CMH|)",
    group_str=group_str_regex,
    # Constraint for focus timepoints - all possible values from the timepoint pairs
    focus_tp="|".join(set([tp for tp_pair in valid_focus_timepoints.values() for tp in tp_pair if tp])),
    # Subject ID constraint - alphanumeric with optional underscores/hyphens
    subject_id="[a-zA-Z0-9_-]+",
    # Explicit treatment/control wildcards (used by regional_contrast rule)
    treatment=f"({'|'.join(sorted({gr['treatment'] for gr in config['analysis']['groups_combinations']}))})",
    control=f"({'|'.join(sorted({gr['control'] for gr in config['analysis']['groups_combinations']}))})",
    
# Function to get sample information from metadata file
def get_sample_info():
    """
    Helper function to retrieve sample information.

    This function is responsible for loading and parsing sample metadata,
    typically from a configuration file or input source. It should return
    information necessary for processing samples in the workflow.

    Returns:
        dict or pandas.DataFrame: Sample information containing metadata
        required for the pipeline execution.
    """
    metadata_path = config["input"]["metadata_path"]
    
    if not metadata_path:
        raise ValueError("metadata_path must be provided in the config file")
    
    # Read metadata file
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    
    # Validate required columns
    required_cols = ["sample_id", "bam_path"]
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        raise ValueError(f"Metadata file is missing required columns: {', '.join(missing_cols)}")
    
    # Create mapping from sample ID to BAM path
    sample_to_bam_map = dict(zip(metadata_df["sample_id"], metadata_df["bam_path"]))
    sample_ids = list(metadata_df["sample_id"])
    
    # Validate that all BAM files exist
    for sample_id, bam_path in sample_to_bam_map.items():
        if not os.path.exists(bam_path):
            raise ValueError(f"BAM file not found for sample {sample_id}: {bam_path}")
    
    return sample_ids, sample_to_bam_map

def _get_mags_by_eligibility(timepoints, groups, eligibility_type):
    """
    INTERNAL: Read the eligibility file for a given timepoint-group combination and return a list of MAG IDs.
    
    WARNING: This function only checks QC eligibility, NOT preprocessing eligibility.
    For rule inputs that need to respect preprocessing eligibility, use get_eligible_mags()
    from dynamic_targets.smk instead.
    
    This function should only be called:
    - Within get_eligible_mags() as a fallback when preprocessing is disabled
    - In preprocessing_eligibility.smk to get initial QC-eligible MAGs
    - In generate_allele_analysis_targets() which runs before preprocessing

    Parameters:
        timepoints (str): The timepoints label
        groups (str): The groups label
        eligibility_type (str or None):
            - "two_sample_unpaired": only return MAG IDs where unpaired_test_eligible is True.
            - "two_sample_paired": only return MAG IDs where paired_test_eligible is True.
            - "lmm": only return MAG IDs where unpaired_test_eligible is True.
            - "all": return MAG IDs that are eligible for any of the tests.
    
    Returns:
        list: MAG IDs that are eligible for the specified test type
    
    Raises:
        FileNotFoundError: If the eligibility file does not exist
    """
    eligibility_file = os.path.join(
        OUTDIR, f"eligibility_table_{timepoints}-{groups}.tsv"
    )
    
    if not os.path.exists(eligibility_file):
        raise FileNotFoundError(
            f"Eligibility file not found: {eligibility_file}. "
            f"Ensure the eligibility_table checkpoint has run for {timepoints}-{groups}."
        )
    
    df = pd.read_csv(eligibility_file, sep="\t")

    if eligibility_type == "two_sample_unpaired" or eligibility_type == "lmm":
        return df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
    elif eligibility_type == "two_sample_paired" or eligibility_type == "cmh":
        return df.loc[df["paired_test_eligible"] == True, "MAG_ID"].tolist()
    # Return MAGs from all eligible columns
    elif eligibility_type == "all":
        # Combine unpaired, paired, and any single-sample eligibility columns.
        single_cols = [
            col for col in df.columns if col.startswith("single_sample_eligible_")
        ]
        return (
            df[
                (df["unpaired_test_eligible"] == True)
                | (df["paired_test_eligible"] == True)
                | (df[single_cols].any(axis=1))
            ]["MAG_ID"]
            .unique()
            .tolist()
        )
    else:
        raise ValueError(
            f"Unknown eligibility type: {eligibility_type}. "
            "Please use 'two_sample_unpaired', 'two_sample_paired', 'cmh', 'lmm' or 'all'."
        )



def _get_single_sample_entries(timepoints, groups):
    """
    INTERNAL: Reads the eligibility file and returns a list of tuples (MAG_ID, group)
    for each column matching 'single_sample_eligible_*' that evaluates to True.
    This allows a MAG to be eligible for multiple single-sample tests.
    
    WARNING: This function only checks QC eligibility, NOT preprocessing eligibility.
    For rule inputs that need to respect preprocessing eligibility, use get_eligible_mags()
    from dynamic_targets.smk with test_type='single_sample' instead.
    """
    eligibility_file = os.path.join(
        OUTDIR, f"eligibility_table_{timepoints}-{groups}.tsv"
    )
    df = pd.read_csv(eligibility_file, sep="\t")
    sample_entries = []
    # Identify all columns with the generic pattern.
    sample_cols = [
        col for col in df.columns if col.startswith("single_sample_eligible_")
    ]

    # If there are no single_sample_eligible columns or all are NaN
    # (which happens for single data_type), return empty list
    if not sample_cols or df[sample_cols].isna().all().all():
        return []

    for _, row in df.iterrows():
        mag = row["MAG_ID"]
        for col in sample_cols:
            if pd.notna(row[col]) and row[col] == True:
                # Extract the sample group (e.g. "control", "fat") from the column name.
                group = col.replace("single_sample_eligible_", "")
                sample_entries.append((mag, group))
    return sample_entries


def get_mags_by_preprocessing_eligibility(timepoints, groups, test_type, group=None):
    """
    Read the preprocessing eligibility file for a given timepoint-group combination and 
    return a list of MAG IDs that have sufficient positions after preprocessing.
    
    This function is used AFTER the preprocessing_eligibility checkpoint runs to determine
    which MAGs should proceed to statistical tests.

    Parameters:
        timepoints (str): The timepoints label (e.g., "pre_post")
        groups (str): The groups label (e.g., "fat_control")
        test_type (str):
            - "two_sample_unpaired": Return MAGs eligible for unpaired two-sample tests
            - "two_sample_paired": Return MAGs eligible for paired two-sample tests
            - "lmm": Return MAGs eligible for LMM analysis
            - "cmh": Return MAGs eligible for CMH analysis
            - "single_sample": Return MAGs eligible for single-sample test (requires group param)
            - "lmm_across_time": Return MAGs eligible for LMM across-time (requires group param)
            - "cmh_across_time": Return MAGs eligible for CMH across-time (requires group param)
        group (str, optional): Required for single_sample, lmm_across_time, cmh_across_time tests
    
    Returns:
        list: MAG IDs that are eligible for the specified test type
    
    Raises:
        FileNotFoundError: If the preprocessing eligibility file does not exist
    """
    # Determine which eligibility file to read based on test type
    # Map test types to the actual eligibility columns:
    #   - two_sample_unpaired, lmm -> two_sample_unpaired_eligible
    #   - two_sample_paired, cmh -> two_sample_paired_eligible
    #   - single_sample, lmm_across_time, cmh_across_time -> single_sample_eligible_{group}
    if test_type in ["two_sample_unpaired", "two_sample_paired", "lmm", "cmh"]:
        eligibility_file = os.path.join(
            OUTDIR, "preprocessing_eligibility", 
            f"preprocessing_eligibility_between_groups_{timepoints}-{groups}.tsv"
        )
        # Map LMM to unpaired eligibility, CMH to paired eligibility
        if test_type in ["two_sample_unpaired", "lmm"]:
            eligible_column = "two_sample_unpaired_eligible"
        else:  # two_sample_paired, cmh
            eligible_column = "two_sample_paired_eligible"
    elif test_type in ["single_sample", "lmm_across_time", "cmh_across_time"]:
        if group is None:
            raise ValueError(f"group parameter is required for test_type '{test_type}'")
        eligibility_file = os.path.join(
            OUTDIR, "preprocessing_eligibility",
            f"preprocessing_eligibility_within_groups_{timepoints}-{groups}.tsv"
        )
        # All within-group tests use single_sample_eligible
        eligible_column = f"single_sample_eligible_{group}"
    else:
        raise ValueError(
            f"Unknown test type: {test_type}. "
            "Use 'two_sample_unpaired', 'two_sample_paired', 'lmm', 'cmh', "
            "'single_sample', 'lmm_across_time', or 'cmh_across_time'."
        )
    
    # Error if file doesn't exist - preprocessing must not have run
    if not os.path.exists(eligibility_file):
        raise FileNotFoundError(
            f"Preprocessing eligibility file not found: {eligibility_file}. "
            f"Ensure preprocessing is enabled and the preprocessing_eligibility checkpoint has run."
        )
    
    df = pd.read_csv(eligibility_file, sep="\t")
    
    if eligible_column not in df.columns:
        logger.warning(
            f"Column '{eligible_column}' not found in {eligibility_file}. "
            "Returning empty list."
        )
        return []
    
    return df.loc[df[eligible_column] == True, "MAG_ID"].tolist()


def parse_metadata_for_timepoint_pairs(timepoints_label, groups_label):
    """
    Parse metadata to identify ancestral-derived sample pairs for dN/dS analysis.
    
    For longitudinal data, this function:
    1. Identifies the derived timepoint from the 'focus' field in config
    2. Treats the other timepoint as ancestral
    3. Matches samples by subject ID across timepoints
    4. Validates exactly 2 samples per subject-timepoint combination
    
    Parameters:
        timepoints_label (str): Timepoint combination label (e.g., "pre_post")
        groups_label (str): Group combination label (e.g., "fat_control")
    
    Returns:
        list: Tuples of (subject_id, ancestral_sample_id, derived_sample_id)
    """
    metadata_path = config["input"]["metadata_path"]
    metadata_df = pd.read_csv(metadata_path, sep="\t")
    
    # Parse timepoint and group information
    if "_" in timepoints_label and DATA_TYPE == "longitudinal":
        timepoint1, timepoint2 = timepoints_label.split("_")
        
        # Determine ancestral and derived timepoints based on focus
        focus_tp = focus_timepoints.get(timepoints_label)
        if not focus_tp:
            raise ValueError(f"No focus timepoint defined for {timepoints_label}")
        
        # The focus is the derived timepoint
        if focus_tp == timepoint1:
            derived_tp, ancestral_tp = timepoint1, timepoint2
        else:
            derived_tp, ancestral_tp = timepoint2, timepoint1
    else:
        raise ValueError(f"dN/dS analysis requires longitudinal data with two timepoints, got: {timepoints_label}")
    
    # Parse groups
    group1, group2 = groups_label.split("_")
    
    # Filter metadata for relevant samples
    filtered_df = metadata_df[
        metadata_df["group"].isin([group1, group2]) &
        metadata_df["time"].isin([ancestral_tp, derived_tp])
    ]
    
    # Group by subject and timepoint to validate sample counts
    subject_timepoint_counts = filtered_df.groupby(["subjectID", "time"]).size()
    
    # Validate exactly 2 samples per subject-timepoint
    for (subject, tp), count in subject_timepoint_counts.items():
        if count != 1:
            raise ValueError(
                f"Expected exactly 1 sample for subject {subject} at timepoint {tp}, "
                f"but found {count} samples"
            )
    
    # Build sample pairs
    sample_pairs = []
    subjects = filtered_df["subjectID"].unique()
    
    for subject in subjects:
        subject_df = filtered_df[filtered_df["subjectID"] == subject]
        
        # Get ancestral and derived samples
        ancestral_samples = subject_df[subject_df["time"] == ancestral_tp]["sample_id"].tolist()
        derived_samples = subject_df[subject_df["time"] == derived_tp]["sample_id"].tolist()
        
        if len(ancestral_samples) == 1 and len(derived_samples) == 1:
            sample_pairs.append((subject, ancestral_samples[0], derived_samples[0]))
        else:
            logger.warning(
                f"Skipping subject {subject}: found {len(ancestral_samples)} ancestral "
                f"and {len(derived_samples)} derived samples"
            )
    
    return sample_pairs


# =============================================================================
# Input Path Helper Functions
# =============================================================================
# These helpers centralize the logic for determining input file paths based on
# data type and configuration options, reducing duplication across rule files.

def get_allele_analysis_input_path(mag_wildcard="{mag}", tp_wildcard="{timepoints}", gr_wildcard="{groups}"):
    """
    Get the appropriate allele analysis input file path based on data type and config.
    
    This helper centralizes the logic for selecting the correct allele frequency file:
    - For single data type: uses filtered or unfiltered single file
    - For longitudinal: uses mean allele frequency changes
    
    Parameters:
        mag_wildcard: MAG ID wildcard string (default: "{mag}")
        tp_wildcard: Timepoints wildcard string (default: "{timepoints}")
        gr_wildcard: Groups wildcard string (default: "{groups}")
    
    Returns:
        str: Path to the appropriate input file
    """
    base_dir = os.path.join(
        OUTDIR,
        "allele_analysis",
        f"allele_analysis_{tp_wildcard}-{gr_wildcard}"
    )
    
    if DATA_TYPE == "single":
        if not config["quality_control"].get("disable_zero_diff_filtering", False):
            return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_no_constant.tsv.gz")
        else:
            return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_single.tsv.gz")
    else:  # longitudinal
        return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_changes_mean.tsv.gz")


def get_allele_freq_cache_path(
    mag_wildcard="{mag}",
    groups_wildcard="{groups}",
    timepoint_wildcard="{timepoint}",
):
    """Path to the per-(MAG, gr_combo, timepoint) allele-frequency Parquet cache file.

    Cache is scoped to gr_combo (not just timepoint) so that the cache rule's
    input can be a single canonical QC file from the same gr_combo — avoiding
    cross-combination QC file dependencies that don't exist at checkpoint time.
    """
    return os.path.join(
        OUTDIR,
        "allele_freq_cache",
        f"{mag_wildcard}_{groups_wildcard}_{timepoint_wildcard}_allele_frequency.parquet",
    )


def get_canonical_qc_file(mag_wildcard, timepoint, groups):
    """Return the single canonical QC TSV path for a (timepoint, gr_combo) cache.

    Selects the first tp_combo in config order that contains ``timepoint`` for
    the given ``groups`` label.  This canonical QC file is the ONLY input to
    ``compute_allele_freq_per_timepoint`` — Snakemake resolves it as a regular
    rule dependency (``qc`` + ``generate_metadata``), so the DAG is valid
    regardless of which combination's eligibility_table checkpoint fired first.

    Why one file is sufficient
    --------------------------
    QC pass/fail for an individual sample depends only on that sample's per-MAG
    breadth and coverage, not on which (tp_combo, gr_combo) triggered the QC
    run.  So the set of passing samples at a given timepoint for a given
    gr_combo is identical across all QC files that cover that (gr_combo,
    timepoint) — any one of them is a correct input for the cache.
    """
    canonical_tp = timepoint_gr_to_canonical_tp.get((timepoint, groups))
    if canonical_tp is None:
        raise ValueError(
            f"No tp_combo found for timepoint={timepoint!r}, groups={groups!r}. "
            "Check that the timepoint and groups appear in the config."
        )
    return os.path.join(
        OUTDIR,
        "QC",
        f"QC_{canonical_tp}-{groups}",
        f"{mag_wildcard}_QC.tsv",
    )




def get_preprocessed_between_groups_path(mag_wildcard="{mag}", tp_wildcard="{timepoints}", gr_wildcard="{groups}"):
    """
    Get the preprocessed between-groups file path based on data type.
    
    Parameters:
        mag_wildcard: MAG ID wildcard string (default: "{mag}")
        tp_wildcard: Timepoints wildcard string (default: "{timepoints}")
        gr_wildcard: Groups wildcard string (default: "{groups}")
    
    Returns:
        str: Path to the preprocessed file
    """
    base_dir = os.path.join(
        OUTDIR,
        "significance_tests",
        f"preprocessed_between_groups_{tp_wildcard}-{gr_wildcard}"
    )
    
    if DATA_TYPE == "single":
        return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_preprocessed.tsv.gz")
    else:  # longitudinal
        return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_changes_mean_preprocessed.tsv.gz")


def get_preprocessed_within_groups_path(mag_wildcard="{mag}", group_wildcard="{group}", 
                                        tp_wildcard="{timepoints}", gr_wildcard="{groups}"):
    """
    Get the preprocessed within-groups file path.
    
    Parameters:
        mag_wildcard: MAG ID wildcard string (default: "{mag}")
        group_wildcard: Group wildcard string (default: "{group}")
        tp_wildcard: Timepoints wildcard string (default: "{timepoints}")
        gr_wildcard: Groups wildcard string (default: "{groups}")
    
    Returns:
        str: Path to the preprocessed within-groups file
    """
    return os.path.join(
        OUTDIR,
        "significance_tests",
        f"preprocessed_within_groups_{tp_wildcard}-{gr_wildcard}",
        f"{mag_wildcard}_{group_wildcard}_allele_frequency_changes_mean_zeros_processed.tsv.gz"
    )

