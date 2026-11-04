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
import pandas as pd
from scipy.stats import chi2
from scipy.stats import f as f_dist
from statsmodels.stats.multitest import multipletests

NOTE_NO_EVENTS = "no_events"
NOTE_NO_EXPOSURE = "no_exposure"
NOTE_NO_RESIDUAL_DF = (
    "no_residual_df"  # tested, but the dispersion could not be estimated
)


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
        The discounted verdict: ``chi2.sf(deviance / dispersion, 1)`` (Test A)
        or ``f.sf(deviance / dispersion, 1, df_resid)`` (Test B); NaN when
        untestable.
    p_value_raw
        The same surprise WITHOUT the cage-agreement discount:
        ``chi2.sf(deviance, 1)``, the textbook Poisson likelihood-ratio test
        that treats every cage as an independent Poisson draw.  Reported
        next to ``p_value`` so a reader sees what requiring the cages to
        agree cost (collaborator request 2026-11-04).  Chi-square for both
        tests: no dispersion is estimated here, so there is nothing for an
        F curve to charge for.  Defined whenever the statistic is (it needs
        no spare cage), so Test B with one cage per group has a raw p but
        no discounted one.
    n
        Cages used (after dropping cages with no exposure in Test A).
    note
        ``""`` when tested, else why not: ``no_events`` / ``no_exposure``;
        ``no_residual_df`` for Test B with one cage per group (no replicate,
        no dispersion estimate: untestable, effect size still reported).
    df_resid
        Residual degrees of freedom the dispersion was estimated on (cages
        minus rates estimated); the F reference's second parameter in Test B.
    """

    rate_ratio: float
    deviance: float
    dispersion: float
    p_value: float
    n: int
    note: str = ""
    df_resid: int = 0
    p_value_raw: float = math.nan


def _untestable(note: str, n: int, rate_ratio: float = math.nan) -> RateTest:
    return RateTest(rate_ratio, math.nan, math.nan, math.nan, n, note, 0)


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
    Such a cage observing an event is impossible by construction (the gene
    tables zero a gene the cage did not read) and raises.

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
    # A cage with expected 0 read none of the gene (or had no events anywhere),
    # so it cannot hold an event of this gene: an observed count there means
    # the gene table and the cage table disagree, which is a bug upstream, not
    # a cage to skip quietly.
    orphan = (e == 0) & (o > 0)
    if orphan.any():
        raise ValueError(
            f"{int(orphan.sum())} cage(s) observe events where expected 0: "
            f"observed {o[orphan].tolist()} with no readable bases in the gene"
        )
    keep = e > 0
    o, e = o[keep], e[keep]
    n = int(len(o))
    if n == 0:
        return _untestable(NOTE_NO_EXPOSURE, 0)
    total_o, total_e = float(o.sum()), float(e.sum())
    if total_o == 0:
        return _untestable(NOTE_NO_EVENTS, n, rate_ratio=0.0)
    # Hand case used throughout: o = [3, 0, 1], e = [0.2, 0.14, 0.5].
    #
    # 1. The rate ratio.  Pool the group's cages: 4 events seen where 0.84
    #    were owed, so the gene collects at 4 / 0.84 = 4.76 times its genome's
    #    pace.  This is the effect size; 1 means "an ordinary gene".
    rr = total_o / total_e
    # 2. The surprise score (likelihood-ratio statistic, "lrt").  Two stories
    #    explain the counts: the boring one, rr = 1 (each cage predicted its
    #    fair share e), and the best-fitting one, rr = 4.76 (each cage
    #    predicted rr * e).  Each story gives the data a probability; the
    #    score is 2 * ln(best / boring).  For Poisson counts that log ratio
    #    collapses to the one line below: the o * ln(o) terms are identical in
    #    both stories and cancel.  Hand case: 2 * (4 ln 4.76 - 3.76 * 0.84)
    #    = 6.17.  Zero means the stories fit equally well; bigger means the
    #    "hot gene" story explains the counts better.
    lrt = 2.0 * (total_o * math.log(rr) - (rr - 1.0) * total_e)
    # 3. The dispersion (cage agreement check).  The score only sees the
    #    group total; it cannot tell one loud cage from several agreeing
    #    cages.  Under the best story each cage is predicted rr * e: 0.95,
    #    0.67, 2.38.  Pearson's X2 adds up (seen - predicted)^2 / predicted
    #    over the cages (5.87) and divides by the spare degrees of freedom,
    #    cages minus the one quantity estimated from them (3 - 1 = 2): 2.94.
    #    About 1 means the cages wobble no more than count luck predicts;
    #    2.94 means one cage shouted (3 where 0.95 was predicted) while
    #    another stayed silent (0 where 0.67 was predicted).  Floored at 1:
    #    cages quieter than luck never make the test MORE confident.
    phi = _dispersion(o, rr * e, n - 1)
    # 4. The verdict.  Discount the surprise by the disagreement (6.17 / 2.94
    #    = 2.10) and ask the chi-square curve with 1 degree of freedom (one
    #    claim on trial: "rr = 1") what fraction of pure-luck scores exceed
    #    it: p = 0.147.  Without the discount the same counts would give
    #    p = 0.013; one loud cage does not carry a gene on its own.
    p = float(chi2.sf(lrt / phi, 1))
    # 5. The raw verdict: the same surprise with no discount, chi2.sf(6.17, 1)
    #    = 0.013 in the hand case.  Reported beside p so the cost of asking
    #    the cages to agree is visible gene by gene.
    raw = float(chi2.sf(lrt, 1))
    return RateTest(rr, lrt, phi, p, n, "", max(n - 1, 0), p_value_raw=raw)


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

    Reference curve: the scaled statistic is referred to F(1, cages - 2),
    the quasi-likelihood F-test, because the dispersion is estimated from
    few cages and the chi-square curve ignores that (see the step-4 comment
    for the simulation numbers).  Test A keeps the chi-square curve: its
    per-cage expected counts are mostly below 1, the dispersion sits at its
    floor, and the F curve there halved the power for no calibration gain.

    Example: four treatment cages with one event each over 100 bases, four
    control cages with none -> r_a 0.01, r_b 0, LRT 2 * 4 ln 2 = 5.55,
    X2 = 0 so dispersion 1, p = F.sf(5.55, 1, 6) = 0.057, rate ratio inf;
    the same pattern over eight cages per group gives p = 0.005.
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
    # Planted case used throughout: four treatment cages with 1 event each
    # over 100 readable bases, four control cages with 0 over 100.
    #
    # 1. The three paces.  Treatment 4 / 400 = 0.010 events per base, control
    #    0 / 400 = 0, and the pooled pace as if the groups were one river,
    #    4 / 800 = 0.005.  Boring story: both groups run at 0.005.  Best
    #    story: treatment 0.010, control 0.
    ra, rb = big_oa / big_la, big_ob / big_lb
    r0 = (big_oa + big_ob) / (big_la + big_lb)

    def term(count: float, rate: float) -> float:
        # Each event is a vote for the best story, worth ln(own pace / pooled
        # pace): an event in a group running twice the pooled pace adds
        # ln 2 = +0.69, one in a group at half the pace adds ln 0.5 = -0.69
        # (logs turn the product of per-event probabilities into a sum, so
        # votes add).  A group with no events casts no votes: 0 * ln(0) is 0
        # in the limit, so return 0 instead of asking for the log of zero.
        return 0.0 if count == 0 else count * math.log(rate / r0)

    # 2. The surprise score: 2 * ln(best / boring) = twice the summed votes.
    #    Planted case: 2 * (4 * ln 2 + 0) = 5.55.  Burst case (treatment
    #    10, 0, 0, 0): 2 * 10 * ln 2 = 13.86, bigger, because the totals are.
    lrt = 2.0 * (term(big_oa, ra) + term(big_ob, rb))
    # 3. The dispersion (cage agreement check).  The score only sees group
    #    totals; this is the one place the cage-as-replicate rule lives.
    #    Under the best story each cage is predicted its group's pace times
    #    its bases: 1 for every treatment cage, 0 for every control cage.
    #    Pearson's X2 = sum of (seen - predicted)^2 / predicted over the
    #    cages, divided by the spare degrees of freedom, cages minus the two
    #    paces estimated (8 - 2 = 6), floored at 1.  Planted case: every cage
    #    sits on its prediction, X2 = 0, dispersion 1 (siblings).  Burst
    #    case: predictions 2.5 each, misses 22.5 + 3 * 2.5 = 30, 30 / 6 = 5
    #    (strangers): one cage had a wild month, the group did not change.
    fitted = np.concatenate([ra * la, rb * lb])
    phi = _dispersion(np.concatenate([oa, ob]), fitted, n - 2)
    # 4. The verdict.  Surprise / dispersion, referred to the F curve with
    #    (1, cages - 2) degrees of freedom: one claim on trial ("the two
    #    paces are equal"), and the dispersion itself was estimated from only
    #    cages - 2 spare cages, which the F curve charges for and the
    #    chi-square curve does not.  Decision 2026-11-03 from a calibration
    #    simulation: with heterogeneous cages the chi-square curve raised
    #    false alarms 12 % of the time at 4 vs 4 cages (9 % at 8 vs 8), the F
    #    curve 5 % (6 %), at a small power cost at the 8-cage floor (0.74 vs
    #    0.80 for a 3x effect).  Planted: 5.55 / 1 on F(1, 6) -> p = 0.057
    #    (chi-square would say 0.0185); the same pattern over 8 vs 8 cages
    #    gives p = 0.005.  Burst: 13.86 / 5 = 2.77 on F(1, 6) -> p = 0.147.
    #    With one cage per group there are no spare cages and no replicate:
    #    untestable (note no_residual_df), never a p-value.
    df_resid = n - 2
    # The raw verdict (no discount, no dispersion estimated, so the plain
    # chi-square curve, never F): planted 5.55 -> 0.0185; burst 13.86 ->
    # 0.0002, which is exactly the number the discount exists to temper.
    raw = float(chi2.sf(lrt, 1))
    if df_resid <= 0:
        # One cage per group: no spare cage to estimate the dispersion from
        # and no replicate at all, so no p-value; the effect size is still
        # worth reporting.
        return RateTest(
            math.inf if rb == 0 else ra / rb,
            lrt,
            math.nan,
            math.nan,
            n,
            NOTE_NO_RESIDUAL_DF,
            0,
            p_value_raw=raw,  # needs no spare cage, so it is still reported
        )
    p = float(f_dist.sf(lrt / phi, 1, df_resid))
    # 5. The effect size: treatment pace over control pace.  inf when control
    #    never had an event, 0 when treatment never did.
    if rb == 0:
        ratio = math.inf
    else:
        ratio = ra / rb
    return RateTest(ratio, lrt, phi, p, n, "", df_resid, p_value_raw=raw)


# ---------------------------------------------------------------------------
# Task 5: the per-gene driver (expected counts, cage floors, BH)
# ---------------------------------------------------------------------------

NOTE_BELOW_FLOOR = "below_cage_floor"


def _group_columns(group: str) -> list[str]:
    return [
        f"n_cages_{group}",
        f"n_de_novo_{group}",
        f"expected_{group}",
        f"rate_ratio_vs_genome_{group}",
        f"dispersion_{group}",
        f"p_value_vs_genome_{group}",
        f"q_value_vs_genome_{group}",
        f"p_value_raw_vs_genome_{group}",
        f"q_value_raw_vs_genome_{group}",
        f"note_{group}",
    ]


def result_columns(groups: tuple[str, str]) -> list[str]:
    """The tests table's columns for a (group a, group b) comparison, in order."""
    a, b = groups
    return (
        ["mag_id", "gene_id", "gene_length", "min_covered_fraction"]
        + _group_columns(a)
        + _group_columns(b)
        + [
            f"rate_{a}",
            f"rate_{b}",
            f"rate_ratio_{a}_vs_{b}",
            "dispersion_between",
            "df_resid_between",
            "p_value_between",
            "q_value_between",
            "p_value_raw_between",
            "q_value_raw_between",
            "note_between",
        ]
    )


def _bh(frame: pd.DataFrame, p_col: str, q_col: str, by_mag: bool) -> None:
    """Benjamini-Hochberg over the tested (non-NaN) p-values of one family, in place.

    One family = one column; genes whose test did not run (NaN p) are left
    out of the correction and keep a NaN q.  ``by_mag`` corrects each genome's
    genes on their own (mirrors ``statistics.fdr_group_by_mag_id``); otherwise
    every gene of the comparison is one family.
    """
    frame[q_col] = np.nan
    keys = frame["mag_id"] if by_mag else pd.Series("all", index=frame.index)
    for _, idx in frame.groupby(keys, sort=False).indices.items():
        rows = frame.index[idx]
        tested = rows[frame.loc[rows, p_col].notna().to_numpy()]
        if len(tested):
            frame.loc[tested, q_col] = multipletests(
                frame.loc[tested, p_col].to_numpy(), method="fdr_bh"
            )[1]


def test_every_gene(
    gene_cage: pd.DataFrame,
    cage_totals: pd.DataFrame,
    groups: tuple[str, str],
    *,
    min_cages: int,
    fdr_by_mag: bool,
) -> pd.DataFrame:
    """Run both tests for every gene of every genome; one row per (genome, gene).

    ``gene_cage`` and ``cage_totals`` are the two tables of
    ``gene_de_novo.gene_tables`` (possibly several genomes stacked).  For a
    gene and a group, the cages that ENTERED (``entered`` True: readable bases
    in the gene above the coverage knob) are the replicates:

    * Test A (``gene_vs_genome``), once per group: observed = the cage's
      events in the gene; expected = the cage's genome-wide events times the
      gene's share of the cage's readable genome, ``n_de_novo_total_c *
      covered_bp_cg / covered_bp_total_c``.  Columns ``*_{group}``.
    * Test B (``group_vs_group``), once: per-cage events and readable bases of
      the gene in group a against group b.  Columns ``*_between`` plus
      ``rate_{a}``, ``rate_{b}`` (pooled events / pooled bases) and
      ``rate_ratio_{a}_vs_{b}``.

    Cage floor (decision 8): Test A for a group needs at least ``min_cages``
    entered cages of that group, Test B needs both groups at the floor;
    below it the note reads ``below_cage_floor`` and the statistics are NaN.
    ``min_covered_fraction`` is the smallest covered fraction among the
    entered cages of both groups (NaN when none entered).

    BH (decision 9): three families, ``p_value_vs_genome_{a}``,
    ``p_value_vs_genome_{b}`` and ``p_value_between``, plus their
    undiscounted twins ``p_value_raw_*`` (collaborator request 2026-11-04:
    the same statistic looked up without the cage-agreement discount, on
    the chi-square curve), each corrected over its tested genes only, per genome when ``fdr_by_mag`` else over the whole
    comparison.

    Worked case (the regression plant's shape, 4 cages per group, a 100-base
    gene in a 1,000-base genome where each treatment cage's only event is in
    that gene): treatment o = [1,1,1,1], e = [0.1]*4 -> rate ratio 10,
    p 8.1e-4 (chi-square reference); control has no events -> ``no_events``;
    between: 4 over 400 vs 0 over 400 -> ratio inf, p 0.057 on F(1, 6)
    (``df_resid_between`` 6); the same pattern over 8 cages per group gives
    p 0.005.
    """
    a, b = (str(g) for g in groups)
    cols = result_columns((a, b))
    gc = gene_cage.copy()
    for col in ("mag_id", "gene_id", "group", "replicate"):
        gc[col] = gc[col].astype(str)
    tot = cage_totals.copy()
    for col in ("mag_id", "group", "replicate"):
        tot[col] = tot[col].astype(str)
    if gc.empty:
        return pd.DataFrame(columns=cols)
    present = set(gc["group"]) | set(tot["group"])
    missing = [g for g in (a, b) if g not in present]
    if missing:
        raise ValueError(
            f"group(s) {missing} not in the gene tables (groups present: {sorted(present)})"
        )
    # Cage totals keyed by (mag, group, replicate); every cage of the gene table
    # must have one, otherwise the expected counts cannot be formed.
    totals = tot.set_index(["mag_id", "group", "replicate"])
    if totals.index.has_duplicates:
        raise ValueError("cage_totals has more than one row for a cage")
    cage_keys = gc.set_index(["mag_id", "group", "replicate"]).index.unique()
    orphan = cage_keys.difference(totals.index)
    if len(orphan):
        raise ValueError(f"cages without a totals row: {list(orphan)[:5]}")

    records = []
    for (mag, gene), sub in gc.groupby(["mag_id", "gene_id"], sort=True):
        rec = {"mag_id": mag, "gene_id": gene}
        rec["gene_length"] = int(sub["gene_length"].iloc[0])
        entered = sub[sub["entered"].astype(bool)]
        rec["min_covered_fraction"] = (
            float(entered["covered_fraction"].min()) if len(entered) else np.nan
        )
        per_group = {}
        for g in (a, b):
            cages = entered[entered["group"] == g]
            o = cages["n_de_novo"].to_numpy(dtype="float64")
            bases = cages["covered_bp"].to_numpy(dtype="float64")
            t = (
                totals.loc[[(mag, g, r) for r in cages["replicate"]]]
                if len(cages)
                else None
            )
            # Fair share: the cage's genome-wide events, times the fraction of
            # its readable genome that this gene is.  A cage whose genome has
            # no readable base at all gets 0 (it was dropped by Task 3 anyway).
            if t is not None:
                denom = t["covered_bp_total"].to_numpy(dtype="float64")
                with np.errstate(divide="ignore", invalid="ignore"):
                    e = np.where(
                        denom > 0,
                        t["n_de_novo_total"].to_numpy(dtype="float64") * bases / denom,
                        0.0,
                    )
            else:
                e = np.zeros(0)
            per_group[g] = (o, bases, e)
            rec[f"n_cages_{g}"] = int(len(cages))
            rec[f"n_de_novo_{g}"] = int(o.sum())
            rec[f"expected_{g}"] = float(e.sum())
            if len(cages) < min_cages:
                rec.update(
                    {
                        f"rate_ratio_vs_genome_{g}": np.nan,
                        f"dispersion_{g}": np.nan,
                        f"p_value_vs_genome_{g}": np.nan,
                        f"p_value_raw_vs_genome_{g}": np.nan,
                        f"note_{g}": NOTE_BELOW_FLOOR,
                    }
                )
            else:
                r = gene_vs_genome(o, e)
                rec.update(
                    {
                        f"rate_ratio_vs_genome_{g}": r.rate_ratio,
                        f"dispersion_{g}": r.dispersion,
                        f"p_value_vs_genome_{g}": r.p_value,
                        f"p_value_raw_vs_genome_{g}": r.p_value_raw,
                        f"note_{g}": r.note,
                    }
                )
        oa, la, _ = per_group[a]
        ob, lb, _ = per_group[b]
        rec[f"rate_{a}"] = float(oa.sum() / la.sum()) if la.sum() > 0 else np.nan
        rec[f"rate_{b}"] = float(ob.sum() / lb.sum()) if lb.sum() > 0 else np.nan
        if len(oa) < min_cages or len(ob) < min_cages:
            rec.update(
                {
                    f"rate_ratio_{a}_vs_{b}": np.nan,
                    "dispersion_between": np.nan,
                    "df_resid_between": 0,
                    "p_value_between": np.nan,
                    "p_value_raw_between": np.nan,
                    "note_between": NOTE_BELOW_FLOOR,
                }
            )
        else:
            r = group_vs_group(oa, la, ob, lb)
            rec.update(
                {
                    f"rate_ratio_{a}_vs_{b}": r.rate_ratio,
                    "dispersion_between": r.dispersion,
                    "df_resid_between": r.df_resid,
                    "p_value_between": r.p_value,
                    "p_value_raw_between": r.p_value_raw,
                    "note_between": r.note,
                }
            )
        records.append(rec)
    out = pd.DataFrame(records)
    # Six BH families: the three discounted ones and their raw twins, each
    # corrected over its own tested genes.
    for g in (a, b):
        _bh(out, f"p_value_vs_genome_{g}", f"q_value_vs_genome_{g}", fdr_by_mag)
        _bh(out, f"p_value_raw_vs_genome_{g}", f"q_value_raw_vs_genome_{g}", fdr_by_mag)
    _bh(out, "p_value_between", "q_value_between", fdr_by_mag)
    _bh(out, "p_value_raw_between", "q_value_raw_between", fdr_by_mag)
    return out[cols].reset_index(drop=True)
