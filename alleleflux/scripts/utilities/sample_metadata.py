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
