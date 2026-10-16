"""MIDAS2's per-site depth filter, applied to ONE sample's profile of ONE MAG.

Why this exists: reads from an organism absent from the catalogue align to the
closest genome present and stack up over the stretches where the two are
similar.  At those positions the "allele frequency" is the mixing ratio of two
organisms, which rises and falls with abundance and does so reproducibly, so it
passes a parallel-change test while meaning nothing about selection.  MIDAS2
removes such positions per sample by comparing a site's depth with the
sample's typical depth for that genome.  This module is that rule, with the
house defaults agreed on 2026-10-15.

The rule (midas2 1.1.3, subcommands/run_snps.py 595-607 and 751;
subcommands/merge_snps.py 386-396), per sample x genome:

    mean_coverage = sum(depth) / count(depth > 0)   with depth = A + C + G + T
    for each site:
        A', C', G', T' = each count if count > min_allele_reads else 0
        depth' = A' + C' + G' + T'
        drop if depth' < site_depth                (checked first)
        drop if depth' / mean_coverage > site_ratio (checked second; '>' not '>=')

Four house differences from MIDAS2:

* ``site_depth`` defaults to 0, not 5.  The goal is to remove positions with
  far MORE reads than the rest of the genome (another genome's reads), not
  thin positions; the floor stays as a parameter for users who want MIDAS2's
  behaviour (``site_depth=5`` reproduces it exactly, proven by brute force in
  the tests).
* ``min_allele_reads`` defaults to 0, not 2: depth' is simply the position's
  coverage (A+C+G+T).  MIDAS2 first zeroes every base with 2 or fewer reads;
  with the floor off that would only nudge depth' by a few reads at the 3x
  boundary, so it is off by default and stays as a parameter.
* The mean coverage is passed IN.  QC computes the covered-sites mean once per
  sample x MAG and stores it as ``average_coverage``; every reader divides by
  that stored number instead of recomputing it from whatever rows it happens
  to have loaded.  ``mean_coverage()`` here is the function QC uses.
* A mean of 0 (no covered site) makes the ratio test vacuous: there is no pile
  to measure against, so nothing is removed by the ceiling.  MIDAS2 never
  reaches this case because its floor already drops every row of such a
  sample; with the floor off the choice has to be made, and "remove nothing"
  is the one that cannot invent a removal.

Dropping is per sample: another sample keeps the site.  Kept rows keep their
original counts by default (``keep_original_counts=True``; the zeroing, if
any, is a throwaway copy that only decides keep/drop).  MIDAS2 instead
carries the zeroed counts downstream, so a kept (18, 2, 0, 0) becomes
(18, 0, 0, 0) with depth 18 when ``min_allele_reads`` is 2; that behaviour is
``keep_original_counts=False`` (``--zero_low_counts``).  ``site_depth_mask``
alone never modifies its input.

Worked example (house defaults, mean_cov = 16 from QC).  Rows (A, C, G, T):
(10,0,0,0): depth' 10, 10 / 16 = 0.6 -> kept.  (40,0,0,0): 40 / 16 = 2.5 ->
kept.  (2,2,2,2): depth' 8, 8 / 16 = 0.5 -> kept unchanged (MIDAS2's own
settings would zero every count, get depth' 0 and drop it as "low").
(60,0,0,0): 60 / 16 = 3.75 -> dropped as "high": this is the pile the filter
exists for.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np
import pandas as pd

BASES = ("A", "C", "G", "T")  # N is never part of the depth (profiles have no N column)

# MIDAS2 1.1.3 defaults, kept here by name so the config template, the CLI help
# and the tests all point at one place.
MIDAS2_SITE_DEPTH = 5
MIDAS2_SITE_RATIO = 3.0
MIDAS2_MIN_ALLELE_READS = 2


@dataclass(frozen=True)
class SiteDepthFilter:
    """The four knobs of the rule.

    enabled          False -> every row is kept (the summary is still computed).
    site_depth       minimum post-zeroing depth; house default 0 = no floor
                     (MIDAS2 DEFAULT_SITE_DEPTH = 5).
    site_ratio       maximum depth' / mean_coverage, '>' strict
                     (MIDAS2 DEFAULT_SITE_RATIO = 3.0).
    min_allele_reads a base count must be GREATER than this to count toward
                     depth'; house default 0 = every read counts (MIDAS2: 2).
    keep_original_counts
                     house modification of MIDAS2, ON by default (decision
                     2026-10-15, as agreed with the collaborator): the zeroing
                     is done on a throwaway copy that only decides keep/drop,
                     and kept rows continue with their original counts.
                     False = MIDAS2 behaviour (the zeroed counts continue
                     downstream); CLI flag ``--zero_low_counts``.

    Example: ``SiteDepthFilter()`` is the house rule;
    ``SiteDepthFilter(site_depth=5, min_allele_reads=2)`` is MIDAS2's
    merge_snps rule verbatim.
    """

    enabled: bool = True
    site_depth: int = 0
    site_ratio: float = MIDAS2_SITE_RATIO
    min_allele_reads: int = 0
    keep_original_counts: bool = True


def _mean_over_covered(depth: np.ndarray) -> float:
    """Mean of the depths that are > 0; 0.0 when none is (MIDAS2's convention)."""
    covered = depth > 0
    return float(depth[covered].mean()) if covered.any() else 0.0


def _depth(counts: pd.DataFrame) -> np.ndarray:
    """Raw per-row depth A+C+G+T as int64 (int64 so a 200k-read pile cannot overflow int32 sums)."""
    return counts[list(BASES)].to_numpy(dtype=np.int64).sum(axis=1)


def mean_coverage(counts: pd.DataFrame) -> float:
    """MIDAS2 mean genome coverage: raw depth averaged over sites with depth > 0.

    ``counts`` needs columns A, C, G, T.  Sites whose four counts sum to 0 are in
    neither the numerator nor the denominator (run_snps.py 604-607).  Returns
    0.0 when no site is covered (run_snps.py 750-751 leaves mean_coverage at 0).
    This is the number QC stores as ``average_coverage``.

    Example: depths 10, 0, 20 -> (10 + 20) / 2 = 15.0.
    """
    depth = _depth(counts)
    covered = depth > 0
    if not covered.any():
        return 0.0
    return float(depth[covered].sum() / covered.sum())


def site_depth_mask(
    counts: pd.DataFrame, flt: SiteDepthFilter, mean_cov: float
) -> tuple[np.ndarray, dict]:
    """Keep-mask (aligned to ``counts.index``) plus a summary, by the MIDAS2 rule.

    Parameters
    ----------
    counts : DataFrame with columns A, C, G, T (other columns are ignored and
        the frame is not modified).
    flt : the four knobs.
    mean_cov : this sample's mean coverage of this MAG over covered sites, as
        QC stored it (``average_coverage``).  Not recomputed here.

    Returns
    -------
    keep : np.ndarray of bool, one per row of ``counts``, True = keep.
    summary : dict with mean_coverage (the value passed in),
        mean_coverage_after_filter (covered-sites mean of the kept rows, over
        the depths they continue with), n_sites_covered (rows with A+C+G+T > 0),
        n_removed_low, n_removed_high, n_kept.  A row
        failing both tests counts once, as "low" (MIDAS2 ``continue``s at the
        first failing test).

    Example (house defaults, mean_cov 6.0): rows (3,0,0,0) x 4 and (18,0,0,0):
    depth' is 3 and 18; 3 / 6 = 0.5 and 18 / 6 = 3.0, neither > 3.0, so all
    five rows are kept.  With ``SiteDepthFilter(site_depth=5)`` the four
    3-read rows fail the floor: keep = [F, F, F, F, T], n_removed_low 4.
    """
    # The mean arrives on trust from a QC table.  NaN would make every
    # comparison False and keep every row silently; a negative number is a
    # corrupt table.  Both must fail loudly here, not downstream.  (Callers
    # should also require the ``average_coverage_including_zeros`` column in
    # the QC table: a table written before it existed stores a different mean
    # under ``average_coverage``.)
    if not np.isfinite(mean_cov) or mean_cov < 0:
        raise ValueError(
            f"site_depth_mask: mean coverage must be a finite, non-negative number, got {mean_cov!r}"
        )
    raw = counts[list(BASES)].to_numpy(dtype=np.int64)
    covered = raw.sum(axis=1) > 0
    n_rows = len(counts)

    if not flt.enabled:
        return np.ones(n_rows, dtype=bool), {
            "mean_coverage": float(mean_cov),
            "mean_coverage_after_filter": _mean_over_covered(raw.sum(axis=1)),
            "n_sites_covered": int(covered.sum()),
            "n_removed_low": 0,
            "n_removed_high": 0,
            "n_kept": int(n_rows),
        }

    # MIDAS2 merge_snps.py 386-390: a base count must be > min_allele_reads to
    # count (with the house default 0 this is the raw depth, i.e. total_coverage).
    # count.  np.where builds a scratch copy; `counts` itself is untouched, so
    # the rows that survive keep their original A, C, G, T.
    depth = np.where(raw > flt.min_allele_reads, raw, 0).sum(axis=1)
    # merge_snps.py 392: the floor, checked first.
    low = depth < flt.site_depth
    # merge_snps.py 396: the ceiling against the sample's mean coverage, checked
    # only on rows that passed the floor.  With mean 0 there is no covered site
    # to compare against, so the ceiling removes nothing (module docstring).
    if mean_cov > 0:
        high = ~low & (depth / float(mean_cov) > flt.site_ratio)
    else:
        high = np.zeros(n_rows, dtype=bool)
    keep = ~(low | high)
    # The covered-sites mean of what continues downstream: the depths the kept
    # rows will carry (zeroed under MIDAS2 settings, original otherwise).  Next
    # to ``mean_coverage`` (the pre-filter reference) it shows how much pile the
    # filter took out of this sample.
    continuing = (
        raw.sum(axis=1)
        if flt.keep_original_counts or flt.min_allele_reads == 0
        else depth
    )
    return keep, {
        "mean_coverage": float(mean_cov),
        "mean_coverage_after_filter": _mean_over_covered(continuing[keep]),
        "n_sites_covered": int(covered.sum()),
        "n_removed_low": int(low.sum()),
        "n_removed_high": int(high.sum()),
        "n_kept": int(keep.sum()),
    }


def apply_site_depth_filter(
    counts: pd.DataFrame, flt: SiteDepthFilter, mean_cov: float
) -> tuple[pd.DataFrame, dict]:
    """Apply the rule to one sample's profile frame, MIDAS2 style.

    Rows that fail the rule are removed.  With the house defaults
    (``min_allele_reads=0``, ``keep_original_counts=True``) rows that pass are
    returned unchanged.  Only with ``keep_original_counts=False`` AND
    ``min_allele_reads > 0`` do kept rows carry the ZEROED base counts (every
    count <= ``min_allele_reads`` set to 0) and, when the frame has a
    ``total_coverage`` column, the zeroed depth: that is MIDAS2's merge_snps
    behaviour.  Every AlleleFlux reader (the cache worker, the ANI loader)
    calls this function.

    Returns a NEW frame (the input is untouched) with a fresh 0..n-1 index and
    the four base columns in their original dtype, plus the summary dict of
    ``site_depth_mask``.  A disabled filter returns the input frame as is.

    Example (``SiteDepthFilter(site_depth=5, min_allele_reads=2,
    keep_original_counts=False)``, mean 16): rows (18, 2, 0, 0) and
    (2, 2, 2, 2) -> the second has depth' 0 < 5 and is dropped; the first is
    kept as (18, 0, 0, 0), total_coverage 18.  With the default
    ``keep_original_counts=True`` it is kept as (18, 2, 0, 0), 20.
    """
    keep, summary = site_depth_mask(counts, flt, mean_cov)
    if not flt.enabled:
        return counts, summary
    out = counts[keep].reset_index(drop=True)
    if flt.min_allele_reads > 0 and not flt.keep_original_counts:
        # MIDAS2 merge_snps.py 386-390 again, this time for keeps: the counts
        # that continue downstream are the zeroed ones.  Cast back so int32
        # cache columns stay int32.
        base = out[list(BASES)]
        zeroed = base.where(base > flt.min_allele_reads, 0)
        for b in BASES:
            out[b] = zeroed[b].astype(base[b].dtype)
        # total_coverage is USED right after this, not just stored: the cache
        # worker's calculate_frequencies divides A, C, G, T by it to make the
        # four *_frequency columns every statistical test reads, and the column
        # is also kept in the parquet for the archival change tables.  So it
        # must follow the zeroed counts: (18, 2, 0, 0) continuing as
        # (18, 0, 0, 0) must divide by 18, not 20.
        if "total_coverage" in out.columns:
            out["total_coverage"] = zeroed.sum(axis=1).astype(
                out["total_coverage"].dtype
            )
    return out, summary


def add_site_depth_filter_args(parser: argparse.ArgumentParser) -> None:
    """Add the four filter flags to a parser; identical on every CLI that reads profiles.

    Defaults are the house defaults (filter on, no floor, ratio 3.0, every
    read counts, original counts kept); ``--site_depth 5 --min_allele_reads 2
    --zero_low_counts`` is MIDAS2's rule.
    The parser is expected to use ``ArgumentDefaultsHelpFormatter`` like the
    other AlleleFlux CLIs, so defaults are not repeated in the help strings.
    """
    parser.add_argument(
        "--disable_site_depth_filter",
        action="store_true",
        help="Do not apply the per-site depth filter to profiles.",
    )
    parser.add_argument(
        "--site_depth",
        type=int,
        default=0,
        help=f"Minimum depth at a site after zeroing low counts; 0 means no floor (MIDAS2 uses {MIDAS2_SITE_DEPTH}).",
    )
    parser.add_argument(
        "--site_ratio",
        type=float,
        default=MIDAS2_SITE_RATIO,
        help="Drop a site when its depth divided by the sample's mean coverage exceeds this.",
    )
    parser.add_argument(
        "--zero_low_counts",
        action="store_true",
        help="MIDAS2 behaviour: at kept positions, base counts not exceeding --min_allele_reads are set to zero for everything downstream. Off = kept positions keep their original counts.",
    )
    parser.add_argument(
        "--min_allele_reads",
        type=int,
        default=0,
        help=f"A base count must exceed this to count toward site depth; 0 means every read counts (MIDAS2 uses {MIDAS2_MIN_ALLELE_READS}).",
    )


def site_depth_filter_from_args(args: argparse.Namespace) -> SiteDepthFilter:
    """Build the dataclass from a namespace parsed with add_site_depth_filter_args."""
    return SiteDepthFilter(
        enabled=not args.disable_site_depth_filter,
        site_depth=args.site_depth,
        site_ratio=args.site_ratio,
        min_allele_reads=args.min_allele_reads,
        keep_original_counts=not args.zero_low_counts,
    )
