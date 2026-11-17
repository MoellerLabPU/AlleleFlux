"""The run metadata sheet's optional per-sample columns and their contract.

Two scripts read the user's metadata sheet directly: ``alleleflux-metadata``
(the first step of every run, seconds of work) and ``alleleflux-baseline-presence``
(after the statistics, hours later).  Both validate the optional columns with
the ONE function here, so a malformed sheet stops the run at the first step
and the two scripts can never disagree about what a date looks like.

Optional columns (spec 2.2, 2.5), blank when the sheet lacks them:

``collection_date``
    The calendar date the sample was collected, ISO 8601 calendar date ONLY:
    ``YYYY-MM-DD`` (e.g. ``2017-01-10``).  The ``cohort`` baseline scope orders
    samples by it.  Any other spelling (``10/01/2017``, ``2017-1-10``,
    ``20170110``, a date with a time) is a malformed sheet and stops the job:
    dates are never guessed, because a day/month mix-up would silently reorder
    the baseline.  An empty cell means "unknown".
``cohort``
    Enrollment batch label, free text, reported as is.

Baseline witnesses (``add_witness_args``, ``cohort_witness_rows``)
------------------------------------------------------------------
Under the ``cohort`` baseline scope ``alleleflux-baseline-presence`` and
``alleleflux-gene-de-novo`` read samples OUTSIDE their comparison as baseline
witnesses.  The workflow names the timepoints and groups it built QC and
allele-frequency caches for (``--witness_timepoints``, ``--witness_groups``);
the two commands share the flags and the row filter defined here, so they
read exactly the same samples and never ask for a cache that was not built.

Pre-treatment timepoints (``checked_pre_treatment_timepoints``)
---------------------------------------------------------------
``--pre_treatment_timepoints`` names the timepoints sampled BEFORE treatment
began.  Under the dated cohort scope a witness at any other timepoint was
sampled after treatment began, and an allele carried by such witnesses only
is labelled "first seen after treatment began" instead of standing variation.
The flag is defined with the witness flags and checked by one function, so
both commands refuse the same mistakes.
"""

from __future__ import annotations

import pandas as pd

META_OPTIONAL = ("collection_date", "cohort")
COLLECTION_DATE_FORMAT = "%Y-%m-%d"  # strptime spelling of YYYY-MM-DD


def validate_collection_dates(meta: pd.DataFrame) -> None:
    """Refuse any ``collection_date`` that is not a real ``YYYY-MM-DD`` date.

    Blank cells (NaN or "") are allowed and mean "date unknown".  Everything
    else must parse with ``COLLECTION_DATE_FORMAT`` exactly: four-digit year,
    two-digit month, two-digit day, hyphens, nothing more.  A non-matching
    value raises ``ValueError`` naming the format and every offending value
    with its sample, so the sheet can be fixed in one go.  A sheet without
    the column passes.

    Example: ``2017-01-10`` passes; ``10/01/2017``, ``2017-1-10``, ``2017-02-30``
    (no such day) and ``2017-01-10 09:30`` all raise.
    """
    if "collection_date" not in meta.columns:
        return
    raw = meta["collection_date"]
    given = raw.notna() & (raw.astype(str).str.strip() != "")
    if not given.any():
        return
    values = raw[given].astype(str).str.strip()
    # pandas accepts some looser spellings under a format (e.g. "2017-1-10");
    # requiring exactly 10 characters pins the zero-padded ISO form.
    parsed = pd.to_datetime(values, format=COLLECTION_DATE_FORMAT, errors="coerce")
    bad = parsed.isna() | (values.str.len() != 10)
    if bad.any():
        offenders = ", ".join(
            f"{s}={d!r}"
            for s, d in zip(meta.loc[values.index[bad], "sample_id"], values[bad])
        )
        raise ValueError(
            f"collection_date must be an ISO calendar date, YYYY-MM-DD (e.g. 2017-01-10), "
            f"or blank; {int(bad.sum())} value(s) are not: {offenders}"
        )


def add_witness_args(parser) -> None:
    """Add the witness flags to a de novo command.

    ``--witness_timepoints`` and ``--witness_groups`` say which samples
    outside the comparison are read; both default to ``None`` = no restriction
    (a run by hand reads every sample of the sheet, as before the flags
    existed).  ``--pre_treatment_timepoints`` says which of the timepoints
    read were sampled before treatment began; default ``None`` = no split.

    One definition for both commands, so the flags mean the same thing in each.
    """
    parser.add_argument(
        "--witness_timepoints",
        nargs="+",
        default=None,
        help=(
            "cohort scope only: a sample OUTSIDE the comparison is read as a "
            "baseline witness only when its timepoint is one of these.  The "
            "workflow passes every timepoint it built allele-frequency caches "
            "at (analysis.baseline_presence.cohort_witnesses decides which).  "
            "Default: any timepoint of the sheet.  The comparison's own samples "
            "are always read."
        ),
    )
    parser.add_argument(
        "--witness_groups",
        nargs="+",
        default=None,
        help=(
            "cohort scope only: a sample OUTSIDE the comparison is read as a "
            "baseline witness only when its group is one of these.  The "
            "workflow passes every group it built QC and caches for.  Default: "
            "any group of the sheet."
        ),
    )
    parser.add_argument(
        "--pre_treatment_timepoints",
        nargs="+",
        default=None,
        help=(
            "cohort scope on a dated sheet only: the timepoints sampled BEFORE "
            "treatment began.  An allele that a pre-treatment baseline witness "
            "carried is standing variation; one that only witnesses at other "
            "timepoints carried is labelled first_seen_after_treatment_began.  "
            "Default: no split (every allele a witness carried is standing "
            "variation).  An error under any other scope, on an undated sheet, "
            "or for a label that is no timepoint of the samples read."
        ),
    )


def checked_pre_treatment_timepoints(
    pre_treatment, baseline_scope: str, sheet_dated: bool, read_timepoints
):
    """Check ``--pre_treatment_timepoints``; return the labels as text, or None.

    The split "standing variation" / "first seen after treatment began" asks
    one question of the baseline witnesses that carry an allele: was any of
    them sampled before treatment began?  That question only exists when
    witnesses can come from several timepoints, which is the cohort scope on
    a fully dated sheet.  Everywhere else the flag could change nothing, so it
    is refused rather than ignored: a silently unused list reads like a split
    that found no on-treatment allele.

    Parameters
    ----------
    pre_treatment : list or None
        The flag's values (``None`` = flag not given).
    baseline_scope : str
        The run's ``--baseline_scope``.
    sheet_dated : bool
        True when every sample of the sheet carries a collection date.
    read_timepoints : iterable
        The timepoints of the samples the command reads (after the witness
        lists were applied).

    Returns
    -------
    ``None`` when the flag was not given (no split); otherwise the labels as
    text, each once, in the order given.

    Raises
    ------
    ValueError
        The flag is given and: the scope is not ``cohort``; or the sheet is
        not fully dated; or a label is not one of ``read_timepoints`` (naming
        the label and the timepoints that are read).

    Example: ``(["t0"], "cohort", True, ["t0", "t1", "tx"])`` -> ``["t0"]``;
    ``(None, "comparison", False, [...])`` -> ``None``;
    ``(["t9"], "cohort", True, ["t0", "t1"])`` -> error naming ``t9``.
    """
    if pre_treatment is None:
        return None
    if baseline_scope != "cohort":
        raise ValueError(
            f"--pre_treatment_timepoints {list(pre_treatment)} only applies under "
            f"--baseline_scope cohort; this run has --baseline_scope {baseline_scope}"
        )
    if not sheet_dated:
        raise ValueError(
            f"--pre_treatment_timepoints {list(pre_treatment)} needs a collection_date "
            "for every sample: without dates the cohort baseline is the comparison's "
            "earlier timepoint only, so no witness is sampled after treatment began"
        )
    # Text on both sides (the sheet is read as text; a label may be typed as a
    # number), each label once, order kept.
    labels = list(dict.fromkeys(str(t).strip() for t in pre_treatment))
    read = [str(t) for t in read_timepoints]
    unknown = [t for t in labels if t not in read]
    if unknown:
        raise ValueError(
            f"--pre_treatment_timepoints: {unknown} is/are not a timepoint of the "
            f"samples read (timepoints read: {sorted(read)})"
        )
    return labels


def cohort_witness_rows(
    meta: pd.DataFrame,
    in_comparison: pd.Series,
    witness_timepoints=None,
    witness_groups=None,
) -> pd.Series:
    """Which rows of the sheet are read under the cohort baseline scope.

    A row is read when it is one of the comparison's own samples
    (``in_comparison``, always kept) or when it may testify: its ``time`` is in
    ``witness_timepoints`` AND its ``group`` is in ``witness_groups``.  A list
    left as ``None`` puts no condition on that column.  Labels are compared as
    text (the sheet is read as text; a list may hold numbers).

    Parameters
    ----------
    meta : DataFrame
        The metadata sheet (or a part of it) with ``time`` and ``group``.
    in_comparison : Series of bool
        True for the comparison's own samples; same index as ``meta``.
    witness_timepoints, witness_groups : iterable or None
        The timepoints and groups the workflow built caches for.

    Returns
    -------
    Series of bool, indexed like ``meta``.

    Example: the comparison is groups A and B at t0 and t1; the sheet also has
    a sample of A at tx and one of C at t0.  With timepoints ``[t0, t1]`` and
    groups ``[A, B]`` neither is read; adding ``tx`` to the timepoints reads
    the first, adding ``C`` to the groups reads the second; with both lists
    ``None`` both are read.
    """
    may_testify = pd.Series(True, index=meta.index)
    if witness_timepoints is not None:
        may_testify &= (
            meta["time"].astype(str).isin({str(t) for t in witness_timepoints})
        )
    if witness_groups is not None:
        may_testify &= meta["group"].astype(str).isin({str(g) for g in witness_groups})
    return (in_comparison.astype(bool) | may_testify).astype(bool)
