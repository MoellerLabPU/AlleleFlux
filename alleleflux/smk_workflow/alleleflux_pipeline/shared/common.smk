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


def get_runtime(rule_name=None):
    """Wall time in minutes (int), with retry scaling. Mirrors :func:`get_time`.

    Used by the native ``snakemake-executor-plugin-slurm`` profile, which expects
    ``runtime`` (integer minutes) — the canonical SLURM-plugin name for what
    cluster-generic calls ``time`` (HH:MM:SS string).

    Rules set BOTH ``time=get_time(...)`` (consumed by cluster-generic via
    ``--time={resources.time}`` in its sbatch template) AND
    ``runtime=get_runtime(...)`` (consumed by plugin-slurm). The unused
    directive is ignored by the active executor; no per-profile rule files
    needed.
    """
    if rule_name:
        base_time = get_resource(rule_name, "time")
    else:
        base_time = config.get("resources", {}).get("time", "24:00:00")

    base_mins = parse_time_to_minutes(base_time)
    retries = _get_retries_for_rule(rule_name)
    step_raw = _get_step_value(rule_name, "time_step", None)

    if retries > 0 and step_raw:
        step_mins = parse_time_to_minutes(step_raw)
        if step_mins > 0:
            return lambda wildcards, attempt: base_mins + (attempt - 1) * step_mins

    return base_mins


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

def check_timepoint_order(meta, earlier, later):
    """Does the sheet's ``collection_date`` agree that ``earlier`` comes before ``later``?

    The config's list order IS the time order for every longitudinal step, so a
    swapped pair silently flips every frequency change.  When the metadata
    sheet carries collection dates (ISO ``YYYY-MM-DD``, the same column the
    cohort baseline scope reads) the order can be checked before any job is
    submitted.  Mice sampled at both timepoints decide: ONE mouse whose
    later-timepoint sample is dated before its earlier one stops the run,
    naming the mice (user ruling 2026-11-02: whether the pair is swapped or
    the sheet's dates are wrong, neither should run).  No subject at both
    timepoints at all raises too: longitudinal AlleleFlux pairs each subject's
    two samples (the allele-frequency step inner-joins on ``subjectID``), so
    such a sheet has nothing to analyse.  No ``collection_date`` column, or
    a blank date among the two timepoints' rows, means no check (``checked``
    False): dates are never guessed.  A date that is not ``YYYY-MM-DD`` raises,
    naming it.

    Returns ``{"checked", "paired"}``.

    Example: m1 pre 2024-01-15 / post 2024-04-15, m2 the same -> checked, 2
    paired; the same sheet asked about ("post", "pre") raises naming both mice;
    a sheet whose pre mice and post mice are different animals raises.
    """
    none = {"checked": False, "paired": 0}
    if "collection_date" not in meta.columns:
        return none
    rows = meta[meta["time"].astype(str).isin([str(earlier), str(later)])].copy()
    dates = rows["collection_date"]
    if dates.isna().any() or (dates.astype(str).str.strip() == "").any():
        return none
    # Strict ISO parse: a day/month mix-up would reorder the check itself.
    parsed = pd.to_datetime(dates.astype(str), format="%Y-%m-%d", errors="coerce")
    if parsed.isna().any():
        bad = sorted(set(dates.astype(str)[parsed.isna()]))
        raise ValueError(
            f"collection_date must be YYYY-MM-DD; cannot read {bad} for timepoints "
            f"{earlier!r} / {later!r}"
        )
    rows["_date"] = parsed.to_numpy()
    rows["time"] = rows["time"].astype(str)
    rows["subjectID"] = rows["subjectID"].astype(str)
    t0 = rows[rows["time"] == str(earlier)].groupby("subjectID")["_date"].min()
    t1 = rows[rows["time"] == str(later)].groupby("subjectID")["_date"].min()
    paired = t0.index.intersection(t1.index)
    if len(paired):
        backward = (t1[paired] < t0[paired]).to_numpy()
        if backward.any():
            bad = sorted(paired[backward].astype(str))
            raise ValueError(
                f"timepoint combination [{earlier}, {later}] lists {later!r} as the later "
                f"timepoint, but {int(backward.sum())} of {len(paired)} mice sampled at "
                f"both have their {later!r} sample dated BEFORE their {earlier!r} sample: "
                f"{bad}.  Either the pair is swapped (exchange the earlier and later "
                "names in the config) or the sheet's collection dates are wrong for "
                "these mice; neither may run"
            )
        return {"checked": True, "paired": int(len(paired))}
    # No subject at both timepoints: the allele-frequency step would inner-join
    # the two timepoints on subjectID and keep nothing, so stop here with the
    # reason instead of hours later with empty tables.
    raise ValueError(
        f"timepoint combination [{earlier}, {later}]: no subject is sampled at both "
        f"timepoints ({len(t0)} subjects at {earlier!r}, {len(t1)} at {later!r}, none "
        "shared).  Longitudinal analysis pairs each subject's two samples; check "
        "that subjectID is spelled the same at both timepoints"
    )


def sheet_fully_dated(meta_sheet):
    """Does every sample of the metadata sheet carry a collection date?

    The ``cohort`` baseline scope orders ALL samples of the sheet by calendar
    date (a witness is a sample dated before the first carrier), so it only
    works on a sheet where every row has a date.  The de novo CLIs make the
    same check on the same sheet and narrow the scope themselves when it
    fails; this is the workflow's copy of that question, asked at parse time
    to decide whether extra timepoints need QC and caches at all.

    Returns True only when the ``collection_date`` column exists and no cell
    is blank.  Whether the dates are well-formed is checked where they are
    used (``check_timepoint_order``, ``extra_cache_timepoints``).

    Example: 4 samples, all with a date -> True; the same sheet with one
    empty cell, or without the column -> False.
    """
    if "collection_date" not in meta_sheet.columns:
        return False
    dates = meta_sheet["collection_date"]
    # The sheet is read with dtype=str, so an empty TSV cell arrives as NaN;
    # a cell holding only spaces arrives as text and is just as blank.
    blank = dates.isna() | (dates.astype(str).str.strip() == "")
    return not bool(blank.any())


def extra_cache_timepoints(meta_sheet, comparison_timepoints):
    """The sheet's timepoints that sit in no comparison, earliest-collected first.

    Under a dated ``cohort`` baseline scope every sample of the sheet dated
    before the first carrier is a baseline witness, whatever its timepoint.  A
    timepoint that is in the sheet but in no ``timepoints_combinations`` entry
    has no QC and no allele-frequency cache unless the workflow asks for them;
    this function names those timepoints so it can.

    Parameters: ``meta_sheet`` is the metadata sheet read as text (columns
    ``sample_id``, ``time``, ``collection_date``); ``comparison_timepoints``
    are the single timepoints of every configured comparison.  Both sides are
    compared as text, so a sheet cell ``5`` and a config value 5 are the same
    timepoint.

    Returns the extra labels ordered by the first collection date of each
    (ties broken by label, so the order never depends on row order).

    Raises ``ValueError`` when: a ``time`` cell is blank (naming the samples);
    an extra label is not a plain token of ASCII letters, digits, ``.`` and
    ``-`` (naming the labels); a date at an extra timepoint is not
    ``YYYY-MM-DD`` (naming the dates).  Only call it on a fully dated sheet
    (``sheet_fully_dated``).

    Example: comparisons use ``a`` and ``b``; the sheet also has ``x`` (first
    sample 2019-06-01) and ``y`` (first sample 2021-01-01) -> ``["x", "y"]``.
    A sheet whose every timepoint is compared -> ``[]``.
    """
    # A blank time cell would become the text "nan" below and turn into a
    # timepoint of its own: stop here, naming the rows.
    blank_time = meta_sheet["time"].isna() | (
        meta_sheet["time"].astype(str).str.strip() == ""
    )
    if blank_time.any():
        raise ValueError(
            "blank time cell in the metadata sheet for sample(s) "
            f"{sorted(meta_sheet.loc[blank_time, 'sample_id'].astype(str))}"
        )
    compared = {str(t) for t in comparison_timepoints}
    times = meta_sheet["time"].astype(str).str.strip()
    extra = sorted(set(times) - compared)
    if not extra:
        return []
    # An extra label becomes a wildcard value and a folder name
    # (QC_{timepoint}, allele_freq_cache/{timepoint}).  "_" is refused because
    # comparison labels are two timepoints joined by "_" and split back on it;
    # anything outside ASCII letters, digits, "." and "-" is refused because it
    # would be a regex operator in the wildcard constraint or an unsafe path
    # segment.  A constraint that silently never matches is worse than a stop.
    bad = [
        t
        for t in extra
        if not (t.isascii() and t.replace(".", "").replace("-", "").isalnum())
    ]
    if bad:
        raise ValueError(
            f"timepoint label(s) {bad} in the metadata sheet cannot be used as a "
            "wildcard value: a timepoint label may hold only ASCII letters, digits, "
            "'.' and '-' (no '_', spaces or slashes).  They sit in no "
            "timepoints_combinations entry but the cohort baseline scope needs QC "
            "and caches there; rename them in the sheet"
        )
    # Strict ISO parse, as in check_timepoint_order (which only read the
    # comparison timepoints' rows): a day/month mix-up would reorder the list.
    rows = times.isin(extra)
    raw = meta_sheet.loc[rows, "collection_date"].astype(str)
    parsed = pd.to_datetime(raw, format="%Y-%m-%d", errors="coerce")
    if parsed.isna().any():
        raise ValueError(
            "collection_date must be YYYY-MM-DD; cannot read "
            f"{sorted(set(raw[parsed.isna()]))} at timepoint(s) {extra}"
        )
    # First collection date per extra timepoint; the label breaks ties.
    first = parsed.groupby(times[rows]).min()
    return sorted(extra, key=lambda t: (first[t], t))


def cohort_caches_every_timepoint(data_type, baseline_scope, use_baseline_presence,
                                  use_gene_de_novo, sheet_dated):
    """Should QC and allele-frequency caches cover EVERY timepoint of the sheet?

    Only when all four hold (user rule 2026-11-13: no extra work unless a de
    novo analysis asks for it):

    * the data are longitudinal (single-timepoint data have no de novo analysis);
    * the baseline scope is ``cohort`` (every other scope reads the
      comparison's own samples only);
    * the sheet is fully dated (on an undated sheet the de novo CLIs narrow
      the scope themselves, so there is nothing extra to read);
    * baseline presence or the gene-level de novo analysis is switched on
      (they are the only readers of the extra caches).

    Example: ``("longitudinal", "cohort", True, False, True)`` -> True;
    the same with ``sheet_dated`` False, scope ``"comparison"``, data type
    ``"single"`` or both analyses off -> False.
    """
    return bool(
        data_type == "longitudinal"
        and baseline_scope == "cohort"
        and sheet_dated
        and (use_baseline_presence or use_gene_de_novo)
    )


# Defined this early (it used to sit with the other de novo helpers further
# down) because the cache-timepoint list below reads its baseline scope before
# the wildcard constraints are declared.  Body unchanged.
def de_novo_rule_params():
    """``analysis.baseline_presence``'s de novo knobs with the house defaults, checked.

    The knobs of the de novo analysis (spec 2.2-2.4), read ONCE so the
    baseline_presence and rising_origin rules can never disagree:

    * ``baseline_scope`` -- what "absent at baseline" is judged against
      (``cohort`` by default: every sample of the run dated before the first
      carrier; needs ``collection_date`` in the metadata);
    * ``cluster_window_bp`` / ``max_de_novo_neighbours`` -- the cluster rule
      (de novo alleles within the window on one contig in one group are one
      block event; more than the tolerated neighbours -> block_event);
    * ``strain_gate`` / ``min_determined_carriers`` -- which strain verdict
      decides whether a de novo allele stands (``carrier`` by default: the
      allele's own carriers' per-mouse turnover verdicts, at least the floor
      of them determined).

    Every gate but ``none`` reads strain-turnover outputs (the per-mouse tables
    for ``carrier``, the classification table for ``group``, both for
    ``both``), so those gates require ``use_strain_turnover`` and the run stops
    at parse time otherwise, not after the statistics.  Unknown words and
    negative numbers stop it too.

    Example: a config with only ``summary`` in the block ->
    ``{"baseline_scope": "cohort", "cluster_window_bp": 1000,
    "max_de_novo_neighbours": 0, "strain_gate": "carrier",
    "min_determined_carriers": 3}``.
    """
    # The two vocabularies, spelled here (not module-level) so the function is
    # self-contained: the slice tests exec it alone, and the CLIs own the words.
    BASELINE_SCOPES = ("individual", "own_group", "comparison", "cohort")
    STRAIN_GATES = ("carrier", "group", "both", "none")
    analysis = config.get("analysis", {}) or {}
    block = analysis.get("baseline_presence", {}) or {}
    p = {
        "baseline_scope": str(block.get("baseline_scope", "cohort")),
        "cluster_window_bp": int(block.get("cluster_window_bp", 1000)),
        "max_de_novo_neighbours": int(block.get("max_de_novo_neighbours", 0)),
        "strain_gate": str(block.get("strain_gate", "carrier")),
        "min_determined_carriers": int(block.get("min_determined_carriers", 3)),
    }
    if p["baseline_scope"] not in BASELINE_SCOPES:
        raise ValueError(
            f"analysis.baseline_presence.baseline_scope must be one of "
            f"{BASELINE_SCOPES}, got {p['baseline_scope']!r}"
        )
    if p["strain_gate"] not in STRAIN_GATES:
        raise ValueError(
            f"analysis.baseline_presence.strain_gate must be one of {STRAIN_GATES}, "
            f"got {p['strain_gate']!r}"
        )
    if p["strain_gate"] != "none" and not analysis.get("use_strain_turnover", False):
        raise ValueError(
            f"analysis.baseline_presence.strain_gate {p['strain_gate']!r} reads the "
            f"strain-turnover outputs (per-mouse tables and/or the classification "
            f"table): set analysis.use_strain_turnover: true (and use_pairwise_ani), "
            f"or strain_gate: none"
        )
    if p["cluster_window_bp"] < 0 or p["max_de_novo_neighbours"] < 0:
        raise ValueError(
            "analysis.baseline_presence.cluster_window_bp and max_de_novo_neighbours "
            "must be >= 0"
        )
    if p["min_determined_carriers"] < 1:
        raise ValueError("analysis.baseline_presence.min_determined_carriers must be >= 1")
    return p


def check_data_type_features(config):
    """Refuse, while the Snakefile is still being read, every setting that cannot
    mean anything for the configured data type.  Returns the checked data type.

    Single-timepoint data has no earlier timepoint, so nothing that compares a
    mouse with its own past can run on it: dN/dS (ancestral vs derived sample),
    regional contrast (mean change over time), baseline presence / de novo
    (absent before, present after) and strain turnover (same-mouse ANI across a
    transition).  Before this guard those features were silently skipped by the
    target generators, so a single-timepoint run with ``use_dnds: True`` finished
    with no dN/dS table and no word about why.  Now the run stops before any job
    is submitted and the message names every offending key.

    Only an EXPLICIT ``True`` is refused.  An absent key is the house default,
    and the target generators already treat absent as off for single data (the
    config template is written for longitudinal data, so its defaults are not
    the user's fault).

    Also checked, for both data types:

    * every ``timepoints_combinations`` entry has the ONE spelling of its data
      type: ``timepoint: "end"`` on single data, ``earlier: "pre"`` + ``later:
      "post"`` on longitudinal data (the later one is the focus; there is no
      focus key).  The list-plus-focus spelling was retired 2026-11-02 because
      it never said which timepoint came first, while every step took list
      order as time order;
    * ``data_type`` itself is one of the two known words (an unknown word used to
      fall through every ``if DATA_TYPE == ...`` branch, including the one that
      picks the output sub-folder);
    * ``use_strain_turnover`` needs ``use_pairwise_ani`` (it reads the ANI table).

    Example: ``{"analysis": {"data_type": "single", "use_dnds": True,
    "timepoints_combinations": [{"timepoint": ["t1"]}]}}`` raises
    ``ValueError`` naming ``use_dnds``; the same config without the ``use_dnds``
    key returns ``"single"``.
    """
    # Vocabulary lives inside the function so the unit tests can slice this
    # function out of the .smk file and exec it alone (see tests/test_data_type_checks.py).
    data_types = ("single", "longitudinal")
    longitudinal_only = (
        "use_dnds",
        "use_regional_contrast",
        "use_baseline_presence",
        "use_strain_turnover",
        "use_gene_de_novo",
    )
    analysis = config["analysis"]
    data_type = analysis.get("data_type")
    if data_type not in data_types:
        raise ValueError(
            f"analysis.data_type must be one of {data_types}, got {data_type!r}"
        )

    # ONE spelling per data type.  Longitudinal: ``earlier`` and ``later``,
    # nothing else (the later one is the focus).  Single: ``timepoint``, one
    # name.  Names may be YAML strings or bare integers (``10mo`` is a string
    # already; ``10`` is read as int and used as the text "10").  Anything
    # else, in particular the retired ``timepoint: [a, b]`` + ``focus`` form,
    # stops here with the expected spelling in the message.
    combos = analysis.get("timepoints_combinations") or []
    if not combos:
        raise ValueError("analysis.timepoints_combinations must list at least one entry")
    want = {"timepoint"} if data_type == "single" else {"earlier", "later"}
    spelling = (
        '- timepoint: "end"'
        if data_type == "single"
        else '- earlier: "pre"\n      later: "post"'
    )
    for tc in combos:
        keys = set(tc) if isinstance(tc, dict) else set()
        named = keys == want and all(
            isinstance(tc[k], (str, int)) and str(tc[k]).strip() for k in want
        )
        if not named or (
            data_type == "longitudinal" and str(tc["earlier"]) == str(tc["later"])
        ):
            raise ValueError(
                f"analysis.data_type is {data_type!r}, so every entry of "
                f"analysis.timepoints_combinations must be spelled\n    {spelling}\n"
                f"(no other keys, no focus key: the later timepoint is the focus; "
                f"two different names for longitudinal data); offending entry: {tc}"
            )

    if data_type == "single":
        # Explicit True only: an absent key is the (longitudinal-written) template
        # default and already counts as off for single data downstream.
        on = [flag for flag in longitudinal_only if analysis.get(flag, False)]
        if on:
            raise ValueError(
                f"analysis.{', analysis.'.join(on)} cannot run on single-timepoint "
                f"data: there is no earlier timepoint to compare against (dN/dS, "
                f"regional contrast, baseline presence / de novo, strain "
                f"turnover and gene de novo all compare a mouse with its own past).  Set the "
                f"flag(s) to false or use data_type: longitudinal."
            )

    # Only ever reached on longitudinal data (on single data use_strain_turnover
    # was already refused above): the turnover table is derived from the ANI
    # table, so asking for one without the other cannot work.
    if analysis.get("use_strain_turnover", False) and not analysis.get(
        "use_pairwise_ani", False
    ):
        raise ValueError(
            "analysis.use_strain_turnover requires analysis.use_pairwise_ani "
            "(it reads the pairwise ANI table)"
        )
    return data_type


def pairwise_ani_settings(config, data_type):
    """Resolve the two pairwise-ANI *mode* knobs with data-type-aware defaults.

    ``analysis.pairwise_ani.pairs`` picks which sample pairs are compared and
    ``store_snp_locations`` which of them also get per-position SNP rows.  Both
    default to ``within_subject`` (same-mouse pairs across timepoints), which is
    the right default for longitudinal data and a meaningless one for
    single-timepoint data: with one sample per mouse there are no same-mouse
    pairs, so ``within_subject`` would compare nothing and ``transitions`` has no
    transitions to match.  On single data the only sensible comparison is every
    mouse against every other (``all``, the between-mouse strain-sharing
    question), so that becomes the default and the same-mouse modes are refused
    rather than quietly producing an empty table.

    Returns ``{"pairs": ..., "store_snp_locations": ...}``; the other keys of the
    block (``min_cov`` etc.) are read elsewhere and ignored here.

    Example: ``{"analysis": {"data_type": "single"}}`` ->
    ``{"pairs": "all", "store_snp_locations": "none"}``;
    ``{"analysis": {"pairwise_ani": {"pairs": "within_subject"}}}`` on single data
    raises ``ValueError``; on longitudinal data it is kept as written.
    """
    block = config["analysis"].get("pairwise_ani", {}) or {}
    if data_type == "single":
        pairs = block.get("pairs", "all")
        store = block.get("store_snp_locations", "none")
        if pairs != "all":
            raise ValueError(
                f"analysis.pairwise_ani.pairs: {pairs!r} needs same-mouse pairs "
                f"across timepoints, but single-timepoint data has one sample per "
                f"mouse; use pairs: all (every pair of samples)"
            )
        if store == "within_subject":
            raise ValueError(
                "analysis.pairwise_ani.store_snp_locations: 'within_subject' "
                "selects same-mouse pairs, of which single-timepoint data has "
                "none; use none or all"
            )
        return {"pairs": pairs, "store_snp_locations": store}
    return {
        "pairs": block.get("pairs", "within_subject"),
        "store_snp_locations": block.get("store_snp_locations", "within_subject"),
    }


# The global data_type used across all Snakemake files, taken from the config
# only after the guard above has refused every setting that cannot work with it.
DATA_TYPE = check_data_type_features(config)
# Pairwise-ANI mode knobs, resolved once; None when the feature is off so a
# config that keeps the template's longitudinal defaults on single data (with
# ANI switched off) is not refused for a block nothing will read.
PAIRWISE_ANI = (
    pairwise_ani_settings(config, DATA_TYPE)
    if config["analysis"].get("use_pairwise_ani", False)
    else None
)
OUTDIR = config["output"]["root_dir"]
# Read defensively so a config with ``use_dnds: False`` (and no ``dnds:`` section)
# still parses.  This value is only consumed when dN/dS targets are generated,
# i.e. when ``use_dnds`` is True — in which case the ``dnds:`` section is expected.
DN_DS_TEST_TYPE = config.get("dnds", {}).get("dn_ds_test_type", "two_sample_paired_tTest")


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

# =============================================================================
# Artifact reuse for permuted / null runs (reuse_from)
# =============================================================================
# A permuted (null) run reuses the real run's *group-independent* expensive
# artifacts — profiles, QC breadth/coverage, and the allele-frequency cache —
# and only recomputes the group-dependent tail (eligibility -> Stage 2 -> stats
# -> scores) on relabeled data.  ``input.reuse_from`` points at the real run's
# data-type output dir (e.g. ``.../alleleflux_output_1/longitudinal``).  When
# set, the cache/QC/metadata path helpers below resolve to that external dir, so
# the profiling/metadata/QC/cache rules drop out of the DAG exactly the way
# USE_EXISTING_PROFILES already drops the profiling rule.  When unset,
# REUSE_DIR == OUTDIR, so non-reuse behaviour is byte-for-byte unchanged.
REUSE_FROM = config["input"].get("reuse_from", "")
USE_REUSE = bool(REUSE_FROM)
if USE_REUSE and not os.path.isdir(REUSE_FROM):
    raise ValueError(
        f"input.reuse_from is set but is not a directory: {REUSE_FROM}"
    )
REUSE_DIR = REUSE_FROM if USE_REUSE else OUTDIR
if USE_REUSE:
    logger.info(f"Reusing profiles/QC/allele-freq cache from: {REUSE_DIR}")


# Permutation mode detection
# ---------------------------
# A permuted (null) run can be in one of three modes, distinguished by config:
#
#   * Orchestrate : permutation.enabled, no ``seed``, no ``permuted_metadata_dir``
#                   -> run_permutations() (workflow.py) fans out one leaf per
#                      seed; the Snakemake pipeline never sees this top config.
#   * Generate    : permutation.seed (an int) is set -> the ``permute_metadata``
#                   rule builds the per-pair sheets under OUTDIR/permuted_metadata
#                   and the group-dependent rules depend on them.
#   * BYO         : permutation.permuted_metadata_dir is set (no ``seed``) ->
#                   read user-supplied sheets directly; the permute_metadata rule
#                   is NOT wired into the DAG.
PERMUTATION = config.get("permutation", {})
PERMUTATION_SEED = PERMUTATION.get("seed")
PERMUTED_METADATA_DIR = PERMUTATION.get("permuted_metadata_dir", "")
GENERATE_PERMUTED = bool(PERMUTATION.get("enabled") and PERMUTATION_SEED is not None)
# Where the per-pair sheets live: generated under OUTDIR, or the BYO directory.
PERMUTED_DIR = (
    os.path.join(OUTDIR, "permuted_metadata")
    if GENERATE_PERMUTED
    else PERMUTED_METADATA_DIR
)

# Null runs switch off the label-blind extras
# --------------------------------------------
# dN/dS, baseline presence, pairwise ANI and strain turnover take their sample
# labels from the real metadata / QC tables and never read a permuted sheet,
# so inside a null run they would mix shuffled significance results with real
# labels: meaningless files and wasted compute.  The statistics (stage 2, the
# tests, CMH, LMM, eligibility) relabel from the permuted sheet and are the
# point of a null run; these four are not.  Forced off here, in the shared
# config layer, so generate-mode leaves AND bring-your-own-sheet runs are both
# covered (decision 2026-10-17).
LABEL_BLIND_EXTRAS = (
    "use_dnds",
    "use_baseline_presence",
    "use_pairwise_ani",
    "use_strain_turnover",
    # The gene-level de novo analysis reads the sample sheet's own labels (its
    # CLI takes no permuted-sheet flag) and its carrier gate reads the strain
    # turnover tables, which this same switch removes from a null run.
    "use_gene_de_novo",
)


def switch_off_label_blind_extras(config, is_permuted_run):
    """Force the five label-blind feature flags to False in a permuted run.

    Mutates ``config["analysis"]`` in place and returns the flags that were ON
    before (for the log line).  A non-permuted run is left untouched.
    ``use_dnds`` defaults to True when absent (see generate_dnds_analysis_targets),
    so an absent key counts as "was on" for that flag only.

    Example: ``{"analysis": {"use_pairwise_ani": True}}`` in a permuted run ->
    returns ``["use_dnds", "use_pairwise_ani"]`` and leaves all five keys False.
    """
    if not is_permuted_run:
        return []
    analysis = config.setdefault("analysis", {})
    was_on = [
        flag for flag in LABEL_BLIND_EXTRAS if analysis.get(flag, flag == "use_dnds")
    ]
    for flag in LABEL_BLIND_EXTRAS:
        analysis[flag] = False
    return was_on


IS_PERMUTED_RUN = GENERATE_PERMUTED or bool(PERMUTED_METADATA_DIR)
_switched_off = switch_off_label_blind_extras(config, IS_PERMUTED_RUN)
if _switched_off:
    logger.info(
        f"Permuted (null) run: switched off {', '.join(_switched_off)} -- these read "
        "real labels, not the permuted sheet, and are not part of the null tail"
    )


def permuted_metadata_path(groups_label):
    """Absolute path to a comparison pair's permuted sheet (generate or BYO).

    Named ``permuted_metadata_<groups_label>.tsv`` (e.g.
    ``permuted_metadata_1D_AL.tsv``) under :data:`PERMUTED_DIR`.
    """
    return os.path.join(PERMUTED_DIR, f"permuted_metadata_{groups_label}.tsv")


def permuted_metadata_flag(groups_label):
    """Return the ``--permuted_metadata`` CLI flag for a group pair, or ``""``.

    Permuted (null) runs relabel the reused cache/QC **in memory** from a
    per-comparison-pair permuted metadata sheet (see
    ``alleleflux-permute-metadata`` and ``relabel_groups_from_metadata``).  In
    generate mode the sheet is produced by the ``permute_metadata`` rule; in BYO
    mode it is the user-supplied file.  Group-dependent rules pass this flag to
    their CLI so each consumer re-derives the permuted labels at load time.
    Returns ``""`` (no flag) when permutation is disabled or no sheet directory
    is configured, so non-permuted runs are unaffected.
    """
    if not PERMUTATION.get("enabled", False):
        return ""
    if not PERMUTED_DIR:
        return ""
    return f"--permuted_metadata {permuted_metadata_path(groups_label)}"


def permuted_metadata_input(groups_label):
    """Return the sheet path as a rule *dependency* — ONLY in generate mode.

    In generate mode the ``permute_metadata`` rule produces the sheet, so every
    consuming rule must declare it as an ``input`` (else Snakemake could
    schedule the consumer before the sheet exists).  In BYO mode the sheet
    already exists on disk and in non-permuted mode there is no sheet, so this
    returns ``[]`` and those runs gain no dependency — behaviour is byte-for-byte
    unchanged.
    """
    return [permuted_metadata_path(groups_label)] if GENERATE_PERMUTED else []

# Profile reuse configuration
# If profiles_path is specified and exists, use existing profiles instead of
# generating new ones.  A reuse_from run defaults profiles to <reuse_from>/profiles
# when profiles_path is not explicitly set (explicit profiles_path still wins).
EXISTING_PROFILES_PATH = config["input"].get("profiles_path", "")
if not EXISTING_PROFILES_PATH and USE_REUSE:
    EXISTING_PROFILES_PATH = os.path.join(REUSE_FROM, "profiles")
USE_EXISTING_PROFILES = bool(EXISTING_PROFILES_PATH and os.path.isdir(EXISTING_PROFILES_PATH))
PROFILES_DIR = EXISTING_PROFILES_PATH if USE_EXISTING_PROFILES else os.path.join(OUTDIR, "profiles")

if USE_EXISTING_PROFILES:
    logger.info(f"Using existing profiles from: {EXISTING_PROFILES_PATH}")
else:
    logger.info(f"Profiles will be generated in: {PROFILES_DIR}")


# =============================================================================
# Sentinel file conventions
# =============================================================================
# Rules that historically used `directory()` outputs are switched to sentinel
# marker files. Snakemake's `directory()` semantics pre-delete the entire
# directory on every re-run, so an interrupted job leaves an empty directory
# which silently invalidates everything downstream (cf. bug where missing
# inputMetadata_{tp} subdirs forced a full pipeline rebuild).
#
# A sentinel file is written ONLY at the very end of a successful job, so:
#   * Interrupted/failed runs leave no sentinel  -> Snakemake re-runs the rule
#   * Successful runs leave the sentinel         -> Snakemake skips the rule
#   * Existing partial outputs are NOT pre-deleted on re-run, so a fresh run
#     can resume cheaply if the script itself is idempotent.
#
# Helpers below build the canonical sentinel path for each rule's output dir.
SENTINEL_PROFILE = ".profile_done"
SENTINEL_METADATA = ".metadata_done"
SENTINEL_QC = ".qc_done"
SENTINEL_DNDS = ".dnds_done"


def profile_sentinel(sample):
    """Path to the sentinel marker for a single-sample profile directory."""
    return os.path.join(PROFILES_DIR, sample, SENTINEL_PROFILE)


def metadata_sentinel(timepoint):
    """Path to the "done" marker of ONE timepoint's inputMetadata folder.

    Metadata is built once per single timepoint (``pre``), not per comparison
    (``pre_post``): a comparison reads the folders of its two timepoints.

    Resolves under REUSE_DIR so a reuse_from run reads the real run's existing
    sentinel (already on disk) and the generate_metadata rule drops from the DAG.

    Example: ``metadata_sentinel("pre")`` ->
    ``<run>/inputMetadata/inputMetadata_pre/.metadata_done``.
    """
    return os.path.join(
        REUSE_DIR, "inputMetadata", f"inputMetadata_{timepoint}", SENTINEL_METADATA
    )


def qc_dir_of(timepoint):
    """The QC folder of ONE timepoint: one ``{mag}_QC.tsv`` per genome of the run.

    The single spelling of the QC layout: the QC marker, the per-genome table
    path and the eligibility rule's ``--qc_dirs`` all build on it.  Under
    REUSE_DIR, like the marker, because QC is reused by a permuted run.

    Example: ``qc_dir_of("pre")`` -> ``<run>/QC/QC_pre``.
    """
    return os.path.join(REUSE_DIR, "QC", f"QC_{timepoint}")


def qc_sentinel(timepoint):
    """Path to the "done" marker of ONE timepoint's QC folder.

    Resolves under REUSE_DIR so a reuse_from run reads the real run's existing
    QC sentinel and the quality_control rule drops from the DAG.

    Example: ``qc_sentinel("pre")`` -> ``<run>/QC/QC_pre/.qc_done``.
    """
    return os.path.join(qc_dir_of(timepoint), SENTINEL_QC)


def comparison_timepoints(timepoints_label):
    """The single timepoints of one comparison label, earlier first.

    Longitudinal comparisons are labelled ``{earlier}_{later}``; single-timepoint
    data use the timepoint itself as the label.  Rules keyed on the comparison
    (eligibility) call this to find the per-timepoint folders they read.

    Example: longitudinal ``"pre_post"`` -> ``["pre", "post"]``; single
    ``"t1"`` -> ``["t1"]``.

    Raises ``ValueError`` when a longitudinal label is not exactly two parts.
    A timepoint label never holds ``_`` (refused at parse time), so that can
    only be a wiring mistake.
    """
    if DATA_TYPE != "longitudinal":
        return [timepoints_label]
    parts = timepoints_label.split("_")
    if len(parts) != 2:
        raise ValueError(
            f"timepoints label {timepoints_label!r} is not EARLIER_LATER "
            f"(split into {parts})"
        )
    return parts


def dnds_sentinel(timepoints, groups, subject_id):
    """Path to the sentinel marker for a per-subject dN/dS output directory."""
    return os.path.join(
        OUTDIR,
        "dnds_analysis",
        f"{timepoints}-{groups}",
        str(subject_id),
        SENTINEL_DNDS,
    )

def transition_labels(combos) -> str:
    """``EARLIER:LATER`` for every two-timepoint combination, space-joined.

    The one spelling of the ``--transitions`` argument that the pairwise-ANI
    rule (``pairs: transitions``) and the strain-turnover rule pass to their
    CLIs, so the two can never disagree.  Single-timepoint entries (no
    ``later``) contribute nothing.  Values are written as given (a numeric
    timepoint name stays ``0:7``).

    Example: ``[{"earlier": "5mo", "later": "10mo"}, {"earlier": "5mo",
    "later": "16mo"}]`` -> ``"5mo:10mo 5mo:16mo"``.
    """
    return " ".join(f"{tc['earlier']}:{tc['later']}" for tc in combos if "later" in tc)


timepoints_labels = []
focus_timepoints = {}

for time_combo in config["analysis"]["timepoints_combinations"]:
    # Shape already checked by check_data_type_features: a single entry carries
    # ``timepoint``, a longitudinal one ``earlier`` + ``later`` (the focus).
    if DATA_TYPE == "single":
        tp = str(time_combo["timepoint"])
        timepoints_labels.append(tp)
        focus_timepoints[tp] = tp
    else:
        earlier, later = str(time_combo["earlier"]), str(time_combo["later"])
        label = f"{earlier}_{later}"
        timepoints_labels.append(label)
        focus_timepoints[label] = later

# The listed order must agree with the sheet's collection dates (when it has
# them): a swapped pair would flip every frequency change silently.  The sheet
# is a required input of every run (the metadata step reads it first), so a
# path that is empty or not on disk stops here rather than one job later.
if DATA_TYPE == "longitudinal":
    _meta_path = config.get("input", {}).get("metadata_path", "")
    if not _meta_path or not os.path.exists(_meta_path):
        raise ValueError(
            f"input.metadata_path must name the sample sheet; got {_meta_path!r}"
            + ("" if _meta_path else " (empty)")
            + ("" if not _meta_path or os.path.exists(_meta_path) else " (no such file)")
        )
    _meta_sheet = pd.read_csv(_meta_path, sep="\t", dtype=str)
    for _label in timepoints_labels:
        _earlier, _later = _label.split("_")
        _order = check_timepoint_order(_meta_sheet, _earlier, _later)
        if not _order["checked"]:
            logger.warning(
                f"timepoint combination {_label}: the metadata sheet has no complete "
                "collection_date column for these timepoints, so the listed order "
                "(earlier first) is taken on trust"
            )
    # _meta_sheet stays alive: the cache-timepoint list below (after
    # unique_timepoints) reads it once more and deletes it there.

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
# The ``allele_freq_cache`` rule writes one Parquet cache file per (MAG, single
# timepoint).  The cache is group-independent and comparison-independent: every
# (timepoints, groups) combination that names the timepoint reuses the same
# file, so each sample's profile is read once however many comparisons use it.
# The cache job reads the QC table of its own timepoint and waits for that
# timepoint's QC marker (qc_sentinel), so it needs no per-comparison lookup.
#
# Example: comparisons a_b, a_c and b_c with three group pairs each would read
# 9 x 2 = 18 (comparison, timepoint) sample sets; the cache reads 3 (a, b, c).

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
# These are the COMPARISON timepoints.  Together with the sheet-only timepoints
# of a cohort-scope run they form CACHE_TIMEPOINTS (below), the values of the
# {timepoint} wildcard: the single-timepoint wildcard of the metadata, QC and
# cache rules, distinct from {timepoints} (the combo label).
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

# CACHE_TIMEPOINTS: every single timepoint that gets QC and an allele-frequency
# cache, and the values the {timepoint} wildcard may take.
#
# Normally that is unique_timepoints (the comparison timepoints).  Under a dated
# cohort baseline scope with a de novo analysis on, the de novo CLIs read EVERY
# sample of the sheet dated before the first carrier as a baseline witness, so
# the sheet's timepoints that sit in no comparison are added (first-collection-
# date order).  No config knob: the sheet is the contract (design D2).
#
# Example: comparisons a_b and a_c; the sheet also holds timepoint x, collected
# before everything else.  Flag off -> ["a", "b", "c"].  Flag on -> ["a", "b",
# "c", "x"], and x gets metadata, QC and caches like the others.
COHORT_CACHES_EVERY_TIMEPOINT = False
EXTRA_CACHE_TIMEPOINTS = []
if DATA_TYPE == "longitudinal":
    _analysis = config.get("analysis", {}) or {}
    _use_baseline_presence = bool(_analysis.get("use_baseline_presence", False))
    _use_gene_de_novo = bool(_analysis.get("use_gene_de_novo", False))
    # The scope is read ONLY when a de novo analysis is on: de_novo_rule_params()
    # also checks the strain gate against use_strain_turnover, and a config with
    # no de novo analysis must never trip that check.
    if _use_baseline_presence or _use_gene_de_novo:
        COHORT_CACHES_EVERY_TIMEPOINT = cohort_caches_every_timepoint(
            DATA_TYPE,
            de_novo_rule_params()["baseline_scope"],
            _use_baseline_presence,
            _use_gene_de_novo,
            sheet_fully_dated(_meta_sheet),
        )
    if COHORT_CACHES_EVERY_TIMEPOINT:
        EXTRA_CACHE_TIMEPOINTS = extra_cache_timepoints(_meta_sheet, unique_timepoints)
    # Say so only when there IS an extra timepoint: every job re-reads this file,
    # so a line that reports "none" would land in every job log for nothing.
    if EXTRA_CACHE_TIMEPOINTS:
        _n_extra_samples = int(
            _meta_sheet["time"].astype(str).str.strip().isin(EXTRA_CACHE_TIMEPOINTS).sum()
        )
        logger.info(
            "cohort scope: QC and caches at every sheet timepoint; extra "
            f"timepoint(s) {EXTRA_CACHE_TIMEPOINTS} ({_n_extra_samples} samples)"
        )
    del _meta_sheet
CACHE_TIMEPOINTS = unique_timepoints + EXTRA_CACHE_TIMEPOINTS

# QC and the allele-frequency cache are both built once per single timepoint
# and are group-independent: the cache of (MAG, timepoint) reads the QC table
# of that same timepoint (get_timepoint_qc_file), and every comparison naming
# the timepoint reuses both.

# =============================================================================
# Wildcard Constraints
# =============================================================================

wildcard_constraints:
    groups=f"({'|'.join(groups_labels)})",
    timepoints=f"({'|'.join(timepoints_labels)})",
    # Every timepoint with QC and a cache: the comparison timepoints, plus the
    # sheet-only ones under the cohort condition (CACHE_TIMEPOINTS above).
    timepoint=f"({'|'.join(CACHE_TIMEPOINTS)})",
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
            - "between_only": MAG IDs eligible for unpaired OR paired tests,
              excluding single-sample (within-group) eligibility.
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
    # Between-group eligibility only — unpaired OR paired, deliberately
    # EXCLUDING single-sample (within-group) columns. Used by allele-analysis
    # target generation when within-group tests are disabled, so that MAGs
    # eligible *only* for within-group tests do not trigger unconsumed
    # allele-analysis / allele-freq-cache jobs.
    elif eligibility_type == "between_only":
        return (
            df[
                (df["unpaired_test_eligible"] == True)
                | (df["paired_test_eligible"] == True)
            ]["MAG_ID"]
            .unique()
            .tolist()
        )
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
            "Please use 'two_sample_unpaired', 'two_sample_paired', 'cmh', "
            "'lmm', 'between_only' or 'all'."
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
    Return MAG IDs eligible for ``test_type`` after preprocessing.

    Uses the canonical Snakemake checkpoint idiom — ``checkpoints.X.get(...).output``
    — to force Snakemake to materialise the eligibility file on disk before
    this function reads it.  This is what makes the function safe to call
    from target-generator functions during DAG planning: if the file does
    not yet exist, the ``.output`` access raises
    ``IncompleteCheckpointException``, which Snakemake catches and re-
    evaluates the DAG after the checkpoint job has run.  A plain
    ``pd.read_csv`` on a manually-constructed path does NOT do this and
    would crash with ``FileNotFoundError`` on cold-start runs whenever the
    output file is eager-listed elsewhere in the DAG.

    Parameters:
        timepoints (str): The timepoints label (e.g., "pre_post")
        groups (str): The groups label (e.g., "fat_control")
        test_type (str): One of
            "two_sample_unpaired" / "two_sample_paired" / "lmm" / "cmh"
              — between-group tests, gated on ``preprocessing_eligibility_between_groups``
            "single_sample" / "lmm_across_time" / "cmh_across_time"
              — within-group tests, gated on ``preprocessing_eligibility_within_groups``
              (require ``group``)
        group (str, optional): Required for the within-group test types.

    Returns:
        list[str]: MAG IDs eligible for the specified test type.
    """
    # Resolve the eligibility file via the canonical Snakemake checkpoint idiom.
    # Accessing ``.output.out_fPath`` on a checkpoint that hasn't completed
    # raises IncompleteCheckpointException, which Snakemake catches and uses
    # to trigger a DAG re-evaluation after the checkpoint job runs.  This is
    # what makes the function safe to call during DAG planning — and why we
    # don't need an explicit os.path.exists guard or a manual path build.
    if test_type in ["two_sample_unpaired", "two_sample_paired", "lmm", "cmh"]:
        # LMM uses the unpaired eligibility column; CMH uses paired.
        eligible_column = (
            "two_sample_unpaired_eligible"
            if test_type in ["two_sample_unpaired", "lmm"]
            else "two_sample_paired_eligible"
        )
        eligibility_file = checkpoints.preprocessing_eligibility_between_groups.get(
            timepoints=timepoints, groups=groups
        ).output.out_fPath
    elif test_type in ["single_sample", "lmm_across_time", "cmh_across_time"]:
        if group is None:
            raise ValueError(f"group parameter is required for test_type '{test_type}'")
        eligible_column = f"single_sample_eligible_{group}"
        eligibility_file = checkpoints.preprocessing_eligibility_within_groups.get(
            timepoints=timepoints, groups=groups
        ).output.out_fPath
    else:
        raise ValueError(
            f"Unknown test type: {test_type}. "
            "Use 'two_sample_unpaired', 'two_sample_paired', 'lmm', 'cmh', "
            "'single_sample', 'lmm_across_time', or 'cmh_across_time'."
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
            # logger.warning(
            #     f"Skipping subject {subject}: found {len(ancestral_samples)} ancestral "
            #     f"and {len(derived_samples)} derived samples"
            # )
            pass
    
    return sample_pairs


# =============================================================================
# Input Path Helper Functions
# =============================================================================
# These helpers centralize the logic for determining input file paths based on
# data type and configuration options, reducing duplication across rule files.

def get_allele_analysis_input_path(mag_wildcard="{mag}", tp_wildcard="{timepoints}", gr_wildcard="{groups}"):
    """
    Get the appropriate allele analysis input file path based on data type and config.

    Returns Parquet paths (1B refactor): all four allele_freq.py outputs now use
    Parquet/Snappy instead of gzip TSV.  Downstream consumers call
    load_allele_freq_inputs() which detects the format by extension.

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
            return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_no_constant.parquet")
        else:
            return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_single.parquet")
    else:  # longitudinal
        return os.path.join(base_dir, f"{mag_wildcard}_allele_frequency_changes_mean.parquet")


def get_allele_freq_cache_path(
    mag_wildcard="{mag}",
    timepoint_wildcard="{timepoint}",
):
    """Path to the per-(MAG, timepoint) allele-frequency Parquet cache file.

    The cache is group-independent (2B refactor): one file per (MAG, timepoint),
    reused across every (timepoint_combination, group_combination) pair that
    includes that timepoint.  The cache job reads the QC table of that same
    timepoint (``get_timepoint_qc_file``).
    """
    return os.path.join(
        REUSE_DIR,
        "allele_freq_cache",
        timepoint_wildcard,
        f"{mag_wildcard}_{timepoint_wildcard}_allele_frequency.parquet",
    )


def get_timepoint_qc_file(mag_wildcard, timepoint):
    """Path of one genome's QC table at ONE timepoint.

    The table lists every sample of that timepoint with its per-sample QC
    verdict; it is header-only when the genome has no profiled sample there.
    The QC step writes one for every genome of the run, so once the QC marker
    of the timepoint exists (``qc_sentinel``) the table exists too.

    Example: ``get_timepoint_qc_file("MAG_7", "pre")`` ->
    ``<run>/QC/QC_pre/MAG_7_QC.tsv``.
    """
    return os.path.join(qc_dir_of(timepoint), f"{mag_wildcard}_QC.tsv")




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



def get_pairwise_ani_output_path(mag_wildcard="{mag}"):
    """Path to the per-MAG pairwise conANI/popANI table (the rule's primary output).

    Group- and timepoint-independent: one file per MAG covering every QC-passing
    sample, because within-subject cross-timepoint comparisons are the point.
    Resolves under OUTDIR (NOT REUSE_DIR): unlike profiles/QC/the allele-freq
    cache, this table is a PRODUCT of the current run, not a reusable upstream
    artifact -- a permuted/null run computes its own or, more sensibly, disables
    the feature via use_pairwise_ani.
    """
    return os.path.join(OUTDIR, "pairwise_ani", f"{mag_wildcard}_pairwise_ani.tsv")


def presence_rule_params():
    """The allele-presence rule, read from config ONCE for every rule that applies it.

    pairwise_ani decides whether two samples share an allele, and
    baseline_presence decides whether an allele was present at the earlier
    timepoint, with the SAME rule: at least ``min_cov`` reads at the position,
    then at least the null-model bar (from ``fdr`` and ``min_base_quality``)
    AND ``min_freq`` of the reads.  Both rules read this one dict, so a changed
    default changes both tables or neither -- they can never drift apart
    silently.  Defaults match the two CLIs' argparse defaults.

    Example: with no ``pairwise_ani`` block in the config ->
    ``{"min_cov": 5, "min_freq": 0.05, "fdr": 1e-6, "min_base_quality": 30}``.
    """
    ani = config["analysis"].get("pairwise_ani", {})
    return {
        "min_cov": ani.get("min_cov", 5),
        "min_freq": ani.get("min_freq", 0.05),
        "fdr": ani.get("fdr", 1e-6),
        # The error model must assume the base-quality floor the profiles were
        # built with -- read from the profiling section, never a separate knob.
        "min_base_quality": config.get("profiling", {}).get("min_base_quality", 30),
    }


# Evaluated once at parse time; rules reference PRESENCE_RULE["min_cov"] etc.
PRESENCE_RULE = presence_rule_params()


def site_depth_filter_params():
    """``quality_control.site_depth_filter`` with the house defaults for every missing key.

    The per-site depth filter (MIDAS2's rule, applied per sample before any
    reader uses a profile): drop a position whose depth exceeds ``site_ratio``
    times the sample's covered-sites mean coverage (QC's ``average_coverage``).
    House defaults (2026-10-15): no floor, ratio 3.0, every read counts toward
    depth, kept positions keep their original counts.  MIDAS2's own values are
    ``site_depth: 5, min_allele_reads: 2, keep_original_counts: false``.

    Example: with no block in the config ->
    ``{"enabled": True, "site_depth": 0, "site_ratio": 3.0,
    "min_allele_reads": 0, "keep_original_counts": True}``.
    """
    block = (config.get("quality_control", {}) or {}).get("site_depth_filter", {}) or {}
    return {
        "enabled": bool(block.get("enabled", True)),
        "site_depth": int(block.get("site_depth", 0)),
        "site_ratio": float(block.get("site_ratio", 3.0)),
        "min_allele_reads": int(block.get("min_allele_reads", 0)),
        "keep_original_counts": bool(block.get("keep_original_counts", True)),
    }


def site_depth_filter_cli_args():
    """The filter flags, identical on every CLI that reads profiles.

    allele_freq_cache, pairwise_ani and baseline_presence all take the flags
    ``add_site_depth_filter_args`` defines; building the string here, once,
    is what keeps the cache and the ANI/baseline tables filtered identically.

    Example: defaults -> ``"--site_depth 0 --site_ratio 3.0 --min_allele_reads 0"``;
    ``enabled: false`` adds ``--disable_site_depth_filter`` in front and
    ``keep_original_counts: false`` adds ``--zero_low_counts`` at the end.
    """
    p = site_depth_filter_params()
    parts = []
    if not p["enabled"]:
        parts.append("--disable_site_depth_filter")
    parts.append(
        f"--site_depth {p['site_depth']} --site_ratio {p['site_ratio']} "
        f"--min_allele_reads {p['min_allele_reads']}"
    )
    if not p["keep_original_counts"]:
        parts.append("--zero_low_counts")
    return " ".join(parts)




def figures_params():
    """The ``analysis.figures`` block with defaults: score figures + run report.

    Example: ``{"analysis": {}}`` -> enabled True, families [] (= every family
    with a score table), tree_path "", format "png", min_mags_per_phylum 3.
    """
    block = config["analysis"].get("figures", {}) or {}
    return {
        "enabled": bool(block.get("enabled", True)),
        "families": list(block.get("families", []) or []),
        "tree_path": block.get("tree_path", "") or "",
        "format": block.get("format", "png"),
        "min_mags_per_phylum": int(block.get("min_mags_per_phylum", 3)),
    }


def groups_for_label(label):
    """``(treatment, control)`` for a groups label: the inverse of ``groups_labels``.

    ``groups_labels`` spells every configured pair as ``f"{treatment}_{control}"``;
    this maps that spelling back.  A lookup, not a split, because a group name
    may itself contain underscores ("high_fat_control" -> ("high_fat", "control")).
    KeyError for a label no configured pair produces (cannot happen for a
    wildcard, which was matched from groups_labels in the first place).
    """
    pairs = {
        f"{gr['treatment']}_{gr['control']}": (gr["treatment"], gr["control"])
        for gr in config["analysis"]["groups_combinations"]
    }
    return pairs[label]


FIGURES = figures_params()  # score heatmaps/boxplots + the run report (analysis.figures)

# Evaluated once at parse time; the three profile-reading rules reference it in params.
SITE_DEPTH_FILTER_ARGS = site_depth_filter_cli_args()


def de_novo_rule_cli_args():
    """The rising-origin flags built from ``de_novo_rule_params``, once.

    Example: defaults -> ``"--cluster_window_bp 1000 --max_de_novo_neighbours 0
    --strain_gate carrier --min_determined_carriers 3"``.  (The baseline scope
    is the baseline_presence rule's own flag and is not in this string.)
    """
    p = de_novo_rule_params()
    return (
        f"--cluster_window_bp {p['cluster_window_bp']} "
        f"--max_de_novo_neighbours {p['max_de_novo_neighbours']} "
        f"--strain_gate {p['strain_gate']} "
        f"--min_determined_carriers {p['min_determined_carriers']}"
    )


# Evaluated once at parse time, only when the feature is on (a config without
# use_baseline_presence must never trip the gate-vs-turnover check).
if config.get("analysis", {}).get("use_baseline_presence", False):
    DE_NOVO_RULE_PARAMS = de_novo_rule_params()
    DE_NOVO_RULE_ARGS = de_novo_rule_cli_args()
else:
    DE_NOVO_RULE_PARAMS = None
    DE_NOVO_RULE_ARGS = ""


def gene_de_novo_params():
    """``analysis.gene_de_novo``'s knobs on top of the shared de novo rules, checked.

    The gene-level analysis (one job per comparison, ``alleleflux-gene-de-novo``)
    asks Phase 2's questions of every rise with the CAGE as the unit.  It
    reuses the shared rules of ``analysis.baseline_presence`` through
    ``de_novo_rule_params`` (baseline scope, cluster window, strain gate) and
    adds three of its own under ``analysis.gene_de_novo``:

    * ``block_rule`` -- ``exclude`` (default: a block of neighbouring de novo
      alleles in one cage is one imported chunk, not counted) or ``once``
      (every block row counts as one event);
    * ``min_gene_covered_fraction`` -- a cage enters a gene only when at least
      this fraction of the gene is readable in it (default 0.0: any readable
      base; the covered fraction is reported either way);
    * ``min_determined_carriers`` -- the carrier gate's floor for THIS
      analysis (default 1: one carrier mouse with a strain verdict decides;
      the site-level analysis keeps its own, default 3).

    A cage has no across-cage strain vote, so ``strain_gate`` ``group`` / ``both``
    stop here naming the block; ``carrier`` still needs ``use_strain_turnover``
    (checked by ``de_novo_rule_params``).

    Example: a config with no ``gene_de_novo`` block and the house
    ``baseline_presence`` defaults -> ``{"baseline_scope": "cohort",
    "cluster_window_bp": 1000, "max_de_novo_neighbours": 0, "strain_gate":
    "carrier", "min_determined_carriers": 1, "block_rule": "exclude",
    "min_gene_covered_fraction": 0.0}``.
    """
    BLOCK_RULES = ("exclude", "once")
    shared = de_novo_rule_params()
    if shared["strain_gate"] in ("group", "both"):
        raise ValueError(
            f"analysis.baseline_presence.strain_gate {shared['strain_gate']!r} is a "
            f"group vote; analysis.gene_de_novo works per cage and accepts only "
            f"carrier or none"
        )
    block = (config.get("analysis", {}) or {}).get("gene_de_novo", {}) or {}
    p = {
        "baseline_scope": shared["baseline_scope"],
        "cluster_window_bp": shared["cluster_window_bp"],
        "max_de_novo_neighbours": shared["max_de_novo_neighbours"],
        "strain_gate": shared["strain_gate"],
        "min_determined_carriers": int(block.get("min_determined_carriers", 1)),
        "block_rule": str(block.get("block_rule", "exclude")),
        "min_gene_covered_fraction": float(block.get("min_gene_covered_fraction", 0.0)),
    }
    if p["block_rule"] not in BLOCK_RULES:
        raise ValueError(
            f"analysis.gene_de_novo.block_rule must be one of {BLOCK_RULES}, "
            f"got {p['block_rule']!r}"
        )
    if not 0.0 <= p["min_gene_covered_fraction"] <= 1.0:
        raise ValueError(
            f"analysis.gene_de_novo.min_gene_covered_fraction must be within [0, 1], "
            f"got {p['min_gene_covered_fraction']}"
        )
    if p["min_determined_carriers"] < 1:
        raise ValueError("analysis.gene_de_novo.min_determined_carriers must be >= 1")
    return p


def gene_de_novo_cli_args():
    """The gene-de-novo rule's knob flags built from ``gene_de_novo_params``, once.

    Example: defaults -> ``"--baseline_scope cohort --cluster_window_bp 1000
    --max_de_novo_neighbours 0 --strain_gate carrier --min_determined_carriers 1
    --block_rule exclude --min_gene_covered_fraction 0.0"``.
    """
    p = gene_de_novo_params()
    return (
        f"--baseline_scope {p['baseline_scope']} "
        f"--cluster_window_bp {p['cluster_window_bp']} "
        f"--max_de_novo_neighbours {p['max_de_novo_neighbours']} "
        f"--strain_gate {p['strain_gate']} "
        f"--min_determined_carriers {p['min_determined_carriers']} "
        f"--block_rule {p['block_rule']} "
        f"--min_gene_covered_fraction {p['min_gene_covered_fraction']}"
    )


# Evaluated once at parse time, only when the feature is on (same guard as
# DE_NOVO_RULE_PARAMS: a config without use_gene_de_novo must never trip the
# gate-vs-turnover check).
# In a permuted (null) run the flag is already False here: it is one of the
# LABEL_BLIND_EXTRAS switched off above, so neither the parameter check (whose
# carrier gate would ask for the strain turnover that null runs lack) nor the
# rule ever runs there.
if config.get("analysis", {}).get("use_gene_de_novo", False):
    GENE_DE_NOVO_PARAMS = gene_de_novo_params()
    GENE_DE_NOVO_ARGS = gene_de_novo_cli_args()
else:
    GENE_DE_NOVO_PARAMS = None
    GENE_DE_NOVO_ARGS = ""


def get_gene_de_novo_tests_path(timepoints="{timepoints}", groups="{groups}"):
    """The one tests table of a comparison: gene_de_novo/{comparison}/gene_de_novo_tests_{comparison}.tsv."""
    return os.path.join(
        OUTDIR,
        "gene_de_novo",
        f"{timepoints}-{groups}",
        f"gene_de_novo_tests_{timepoints}-{groups}.tsv",
    )


def get_strain_turnover_output_path(mag_wildcard="{mag}"):
    """Per-MAG strain-turnover table (the strain_turnover rule's primary output); OUTDIR, like ANI."""
    return os.path.join(OUTDIR, "strain_turnover", f"{mag_wildcard}_strain_turnover.tsv")


def get_replacement_classification_path():
    """The one all-MAG classification table the enrichment filter reads."""
    return os.path.join(OUTDIR, "strain_turnover", "replacement_classification.tsv")


# The baseline-presence summary family, read once: it names the input summary
# file, the source directory and the output stem.
BASELINE_PRESENCE_FAMILY = config["analysis"].get("baseline_presence", {}).get("summary", "two_sample_paired")
# Which test's rows of that summary are annotated, also read once.  Default: the
# family's t-test ("two_sample_paired" -> "two_sample_paired_tTest").  The rule's
# --test_type and the output stem below BOTH use this constant, so a config
# without a baseline_presence block gives them the same value instead of one of
# them failing with a KeyError while the Snakefile is being read.
BASELINE_PRESENCE_TEST_TYPE = config["analysis"].get("baseline_presence", {}).get(
    "test_type", f"{BASELINE_PRESENCE_FAMILY}_tTest"
)


def get_baseline_presence_stem(timepoints="{timepoints}", groups="{groups}"):
    """Output stem for one comparison: {comparison}_{family}_{statistic}_baseline_presence.

    Mirrors ``baseline_presence.output_label``: the family prefix is stripped
    from the test_type so "two_sample_paired" + "two_sample_paired_tTest" gives
    "two_sample_paired_tTest", while "lmm" + "LMM_abs" gives "lmm_LMM_abs".
    """
    family = BASELINE_PRESENCE_FAMILY
    test_type = BASELINE_PRESENCE_TEST_TYPE
    prefix = f"{family}_"
    stat = test_type[len(prefix):] if test_type.startswith(prefix) else test_type
    return os.path.join(OUTDIR, "baseline_presence", f"{timepoints}-{groups}_{family}_{stat}_baseline_presence")


def get_all_qc_files_for_mag(mag_wildcard="{mag}"):
    """One genome's QC table at every COMPARISON timepoint -- the ANI sample gate.

    Unlike the allele-frequency cache (which wants ONE timepoint's QC),
    pairwise ANI wants every sample judged usable at any timepoint of any
    comparison; the CLI unions them.  One table per entry of
    ``unique_timepoints``: the same samples the per-comparison tables used to
    hold between them, each sample now in exactly one table.  Timepoints that
    sit in no comparison (the extra cache timepoints of a cohort-scope run) are
    left out, as before: ANI compares the samples of the comparisons.

    Example: comparisons ``pre_post`` and ``pre_end`` ->
    ``[<run>/QC/QC_pre/MAG_7_QC.tsv, .../QC_post/..., .../QC_end/...]``.
    """
    return [get_timepoint_qc_file(mag_wildcard, tp) for tp in unique_timepoints]
