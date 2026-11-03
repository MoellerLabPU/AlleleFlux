"""Two Poisson rate tests for de novo events per gene, cage as the replicate.

The question (plan decision 1).  Each cage reports how many new alleles
landed in a gene and how many bases of that gene it could search.  Counts
of 0, 1, 2 over a known amount of searched sequence are Poisson material,
so both tests are likelihood-ratio tests of Poisson counts with an
exposure offset, in closed form (one or two rate parameters, so no fit can
fail), with one honest concession to cages that disagree with each other:
an overdispersion scale.

* ``gene_vs_genome`` (Test A, within one group): is this gene a hotspot of
  its own genome?  Per cage, the expected count is the cage's genome-wide
  event rate times the gene's readable bases; one rate ratio across the
  group's cages; the test asks whether that ratio is 1.
* ``group_vs_group`` (Test B, across groups): does the gene collect events
  per readable base faster in one group than in the other?  One rate per
  group, pooled under the null.

Overdispersion.  A Poisson model trusts its counts exactly; a cage with a
burst of events would make it overconfident.  Pearson's X2 over the
residual degrees of freedom estimates how much noisier the cages are than
Poisson predicts (the quasi-Poisson idea); the LRT statistic is divided by
that scale, floored at 1 (cages quieter than Poisson never make a test
MORE confident), and read off chi-square with 1 degree of freedom.  A
random cage effect was deliberately not fitted: with four cages per group
its variance estimate is unusable (the LMM lesson of 2026-10).

Worked case (Test A): three cages observe 3, 0, 1 events where 0.2, 0.14
and 0.5 were expected.  Rate ratio 4 / 0.84 = 4.76; LRT 6.17; the cages
disagree (X2 5.87 over 2 df -> scale 2.94), so the scaled statistic is 2.10
and p = 0.147.  Without the scale the same numbers would give p = 0.013:
one cage with 3 events where 0.2 were expected is less convincing when
another cage with a similar expectation saw none.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.stats import chi2

NOTE_NO_EVENTS = "no_events"
NOTE_NO_EXPOSURE = "no_exposure"


@dataclass(frozen=True)
class RateTest:
    """One test's result.

    rate_ratio
        Test A: the gene's event rate over its genome's (``sum(o) / sum(e)``);
        Test B: group a's rate over group b's (``inf`` when b has no events,
        ``0`` when a has none).
    deviance
        The likelihood-ratio statistic, 2 * (log-likelihood ratio), before
        scaling.
    dispersion
        ``max(1, Pearson X2 / df_resid)``; 1 when there is no residual df.
    p_value
        ``chi2.sf(deviance / dispersion, 1)``; NaN when untestable.
    n
        Cages used (after dropping cages with no exposure in Test A).
    note
        ``""`` when tested, else why not: ``no_events`` / ``no_exposure``.
    """

    rate_ratio: float
    deviance: float
    dispersion: float
    p_value: float
    n: int
    note: str = ""


def _untestable(note: str, n: int, rate_ratio: float = math.nan) -> RateTest:
    return RateTest(rate_ratio, math.nan, math.nan, math.nan, n, note)


def _as_float(name: str, values) -> np.ndarray:
    arr = np.asarray(values, dtype="float64").ravel()
    if np.isnan(arr).any() or (arr < 0).any():
        raise ValueError(f"{name} must be non-negative numbers, got {values}")
    return arr


def _dispersion(observed: np.ndarray, fitted: np.ndarray, df_resid: int) -> float:
    """Pearson X2 / df_resid, floored at 1; cells with fitted 0 contribute 0."""
    if df_resid <= 0:
        return 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        cells = np.where(fitted > 0, (observed - fitted) ** 2 / fitted, 0.0)
    return max(1.0, float(cells.sum()) / df_resid)


def gene_vs_genome(observed, expected) -> RateTest:
    """Test A: observed events per cage against the cage's expected count.

    ``expected[c]`` is ``n_de_novo_total_c * covered_bp_gene_c /
    covered_bp_total_c`` (the gene's fair share of the cage's events).  The
    model is Poisson with log(expected) as the offset and one free
    parameter, the rate ratio; its maximum-likelihood value is
    ``rr = sum(o) / sum(e)`` and the likelihood-ratio statistic against
    ``rr = 1`` is ``2 * (sum(o) * ln rr - (rr - 1) * sum(e))``.  Cages with
    ``expected == 0`` (nothing readable in the gene, or a genome with no
    events) are dropped first: they carry no information about this gene.

    Untestable: no events at all (``rate_ratio 0``, note ``no_events``; a
    gene nobody hit is not "p = 1", it is unmeasured) or no exposure
    (note ``no_exposure``).

    Example: o = [3, 0, 1], e = [0.2, 0.14, 0.5] -> rr 4.76, LRT 6.17,
    dispersion 2.94, p 0.147 (see the module docstring).
    """
    o = _as_float("observed", observed)
    e = _as_float("expected", expected)
    if len(o) != len(e):
        raise ValueError(
            f"observed and expected differ in length: {len(o)} vs {len(e)}"
        )
    keep = e > 0
    o, e = o[keep], e[keep]
    n = int(len(o))
    if n == 0:
        return _untestable(NOTE_NO_EXPOSURE, 0)
    total_o, total_e = float(o.sum()), float(e.sum())
    if total_o == 0:
        return _untestable(NOTE_NO_EVENTS, n, rate_ratio=0.0)
    rr = total_o / total_e
    # Null model: rr = 1 (fitted = e); alternative: fitted = rr * e.  The
    # Poisson log-likelihood difference reduces to this one line because the
    # o * ln(o) terms cancel between the two models.
    lrt = 2.0 * (total_o * math.log(rr) - (rr - 1.0) * total_e)
    phi = _dispersion(o, rr * e, n - 1)
    p = float(chi2.sf(lrt / phi, 1))
    return RateTest(rr, lrt, phi, p, n)


def group_vs_group(observed_a, exposure_a, observed_b, exposure_b) -> RateTest:
    """Test B: one event rate per group, cages as the replicates.

    ``observed_*`` are the per-cage event counts in the gene and
    ``exposure_*`` the per-cage readable bases of the gene.  Rates
    ``r_a = O_a / L_a``, ``r_b = O_b / L_b`` (group totals), pooled
    ``r_0 = (O_a + O_b) / (L_a + L_b)``; the likelihood-ratio statistic for
    "one rate" against "two rates" is ``2 * (O_a ln(r_a / r_0) + O_b
    ln(r_b / r_0))``, a zero count contributing 0.  Pearson X2 sums over
    every cage against its group's fitted ``r_g * L_cg``; df_resid is the
    number of cages minus 2.

    ``rate_ratio`` is ``r_a / r_b``: ``inf`` when only a has events, ``0``
    when only b has.  Untestable: no events in either group (``no_events``)
    or a group with no readable bases at all (``no_exposure``).

    Example: four treatment cages with one event each over 100 bases, four
    control cages with none -> r_a 0.01, r_b 0, LRT 2 * 4 ln 2 = 5.55,
    X2 = 0 so dispersion 1, p = 0.0185, rate ratio inf.
    """
    oa, la = _as_float("observed_a", observed_a), _as_float("exposure_a", exposure_a)
    ob, lb = _as_float("observed_b", observed_b), _as_float("exposure_b", exposure_b)
    if len(oa) != len(la) or len(ob) != len(lb):
        raise ValueError(
            f"observed and exposure differ in length: a {len(oa)} vs {len(la)}, "
            f"b {len(ob)} vs {len(lb)}"
        )
    n = int(len(oa) + len(ob))
    big_oa, big_la, big_ob, big_lb = (float(v.sum()) for v in (oa, la, ob, lb))
    if big_la == 0 or big_lb == 0:
        return _untestable(NOTE_NO_EXPOSURE, n)
    if big_oa + big_ob == 0:
        return _untestable(NOTE_NO_EVENTS, n)
    ra, rb = big_oa / big_la, big_ob / big_lb
    r0 = (big_oa + big_ob) / (big_la + big_lb)

    def term(count: float, rate: float) -> float:
        # 0 * ln(0 / r0) is 0 in the limit: a group with no events adds nothing.
        return 0.0 if count == 0 else count * math.log(rate / r0)

    lrt = 2.0 * (term(big_oa, ra) + term(big_ob, rb))
    fitted = np.concatenate([ra * la, rb * lb])
    phi = _dispersion(np.concatenate([oa, ob]), fitted, n - 2)
    p = float(chi2.sf(lrt / phi, 1))
    if rb == 0:
        ratio = math.inf
    else:
        ratio = ra / rb
    return RateTest(ratio, lrt, phi, p, n)
