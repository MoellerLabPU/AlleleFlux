#!/usr/bin/env python3
"""Standing variation or de novo, for the alleles that ROSE  (``alleleflux-rising-origin``).

Reads the ``*_baseline_presence_rising.tsv`` tables written by
``alleleflux-baseline-presence`` (one row per significant site x allele x
group the allele rose in) and counts, per comparison x group, how many rising
alleles were already present at the earlier timepoint and how many were not.

Origin categories (``origin``), from the rising table's ``origin_any_mouse``
with the timepoint names taken out so every comparison shares one vocabulary:

    standing_variation       present in >= 1 covered earlier sample ("any" rule)
    de_novo                  covered earlier samples, the allele in none of them
    de_novo_below_detection  only a few reads at the earlier timepoint, under
                             the detection bar
    new_strain               de novo or de novo below detection, but the strain
                             gate says the allele most likely arrived with a new
                             strain, not by mutation: the group's across-cage
                             vote says REPLACED (gate ``group``), or more than
                             half of the allele's own carriers changed strain
                             (gate ``carrier``)
    strain_undetermined      de novo or de novo below detection, but the strain
                             gate could not say: too few carriers with a
                             per-mouse verdict, or a group vote that is neither
                             replaced nor not_replaced (too few voters, no
                             voters, not in the table, a tie).  Never counted
                             as de novo (user ruling 2026-10-23)
    block_event              de novo (either flavour) but with more than
                             ``--max_de_novo_neighbours`` other de novo alleles
                             within ``--cluster_window_bp`` on the same contig,
                             same group, any age: one imported or recombined
                             chunk, never counted as independent mutations
                             (spec 2.3; ``block_id`` names the block)
    not_present_later        no covered later sample of the group shows it
    baseline_not_covered     no earlier sample deep enough to judge

The first four are verdicts (``CALLED``); the last four are not.  The
cluster rule runs BEFORE any strain gate: whether alleles sit together is a
fact about positions, so it is judged on every candidate before a gate can
relabel some of them (a gate that ran first would hide neighbours and could
leave one member of a block looking like a lone de novo allele).

Strain verdicts (``--classification``, optional)
------------------------------------------------
One or more ``replacement_classification.tsv`` tables (alleleflux-replacement-
classification), one per voter floor (read from ``min_voters``).  Every rising
allele is joined to its MAG's verdict for that period, once per floor (floors
are never mixed):

* within-group rows (``group_analyzed`` set) take the verdict of the group the
  allele rose in -- the control arm of a comparison has its own verdict;
* between-group rows take the dirtier of the comparison's two arms
  (``combine_arms``): one replaced arm is enough to explain a difference.

A MAG x group x period with no row in a table was never classified: its status
is ``not_in_strain_table`` and its class ``unverified``, never a guess.
``attach_strain`` only REPORTS (``strain_status``, ``strain_class``); which
verdict changes an origin is the gate's business.

Carrier verdicts (``--turnover_dir``, optional; spec 2.4)
--------------------------------------------------------
The per-mouse ``{mag}_strain_turnover.tsv`` tables of alleleflux-strain-turnover
(one row per mouse and transition, ``dominant_strain_change`` True / False /
blank).  For every rising allele the mice carrying it at the later timepoint
(``{later}_mice_allele_present`` from the rising table) are looked up for the
allele's period: ``carrier_strain_status`` is ``new_strain_carriers`` when more
than half of the determined carriers changed, ``same_strain_carriers``
otherwise, ``strain_undetermined`` when fewer than ``--min_determined_carriers``
could be judged (or the MAG has no turnover table).  The thresholds live in the
turnover step; nothing is recomputed here.

The gate (``--strain_gate``; ``apply_strain_gate``, the ONLY place an origin
changes after the cluster rule) relabels the two de novo flavours and nothing
else:

* ``carrier`` (default): from ``carrier_strain_status``;
* ``group``: from ``strain_class``: replaced -> new_strain, not_replaced -> kept,
  anything else -> strain_undetermined;
* ``both``: an allele must pass both; one undetermined verdict makes it
  ``strain_undetermined`` (not explained), otherwise one ``replaced`` /
  ``new_strain_carriers`` makes it ``new_strain`` (strain_undetermined >
  new_strain > kept);
* ``none``: no relabel.

Whatever the gate, every row carries BOTH verdicts (``strain_status`` /
``strain_class`` and ``carrier_strain_status`` + counts) whenever their inputs
were given; the unused one is reporting only, there to be compared, never read
by the gate.

Outputs (``--output_dir``)
--------------------------
* ``rising_origin_alleles.tsv`` -- every rising allele with its ``origin``,
  ``block_id`` / ``n_de_novo_neighbours``, and, when the inputs were given,
  ``floor`` / ``strain_status`` / ``strain_class`` (group verdict) and
  ``carrier_strain_status`` / ``n_carriers`` / ``n_carriers_determined`` /
  ``n_carriers_strain_changed`` (carrier verdict).
* ``rising_origin_summary.tsv`` -- one row per floor x comparison x group x
  strain class x origin, zero counts included: ``n_alleles``, ``n_mags``,
  ``n_rising`` / ``n_called`` (the block's totals), ``pct_of_all`` (of every
  rising allele in the block) and ``pct_of_called`` (of the four verdicts; NaN
  for the two non-verdicts).  Strain classes: each class present plus ``all``;
  without verdicts only ``all`` and no floor.  Comparisons and groups are
  never pooled: the control appears once per comparison.

Example (group 40, 5 rising alleles: standing 1, de novo 1, de novo below
detection 1, not present later 1, baseline not covered 1; no verdicts):
standing_variation n=1, pct_of_all 20.0, pct_of_called 33.3.
"""

import argparse
import glob
import logging
import os
import re

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.baseline_presence import parse_comparison
from alleleflux.scripts.analysis.de_novo_rules import (
    BLOCK_EVENT,
    CARRIER_NEW_STRAIN,
    CARRIER_SAME_STRAIN,
    CARRIER_UNDETERMINED,
    assign_blocks,
    carrier_strain_verdicts,
)
from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

# Origin categories, in report order; the first four are verdicts.
ORIGINS = [
    "standing_variation",
    "de_novo",
    "de_novo_below_detection",
    "new_strain",
    BLOCK_EVENT,  # not a verdict: a cluster of de novo alleles counted as one event
    CARRIER_UNDETERMINED,  # not a verdict: the strain gate could not say
    "not_present_later",
    "baseline_not_covered",
]
CALLED = ORIGINS[:4]
DE_NOVO = ("de_novo", "de_novo_below_detection")  # the only origins a gate relabels
STRAIN_GATES = ("carrier", "group", "both", "none")
# Precedence of a gate verdict under ``both`` (user ruling 2026-10-26): an
# allele is only EXPLAINED (new_strain) when both gates could speak and at
# least one found a swap; if either gate is undetermined the allele is removed
# from the de novo count as undetermined, and the other gate's verdict stays
# readable in its own column.  None = the gate keeps the origin.
_GATE_RANK = {None: 0, "new_strain": 1, CARRIER_UNDETERMINED: 2}

# Raw verdicts a classification table may carry, and the class each maps to.
# Anything else is a contract violation and raises (never guessed).
STATUS_TO_CLASS = {
    "not_replaced": "not_replaced",
    "replaced": "replaced",
    "tie": "tie",  # only when the table was built with tie = unresolved
    "too_few_voters": "unverified",
    "no_voters": "unverified",
    "not_in_strain_table": "unverified",  # this script's label for a missing row
}

# Columns read from a rising table.
REQUIRED = [
    "mag_id",
    "contig",
    "position",
    "gene_id",
    "group_analyzed",
    "allele",
    "rose_in",
    "frequency_change",
    "n_replicates",
    "n_mice",
    "origin_any_mouse",
]
BLOCK = ["floor", "comparison", "period", "age", "test_type", "rose_in", "strain_class"]
SUMMARY_COLUMNS = BLOCK + [
    "origin",
    "n_alleles",
    "n_mags",
    "n_rising",
    "n_called",
    "pct_of_all",
    "pct_of_called",
]

# {tp1}_{tp2}-{g1}_{g2}_{test label}_baseline_presence_rising.tsv
_NAME = re.compile(
    r"^(?P<comparison>[^_-]+_[^_-]+-[^_]+_[^_]+)_(?P<test>.+)_baseline_presence_rising\.tsv$"
)


def parse_rising_name(path: str) -> tuple[str, str]:
    """(comparison, test label) from a rising table's file name.

    Example: ``5mo_16mo-40_AL_single_sample_tTest_baseline_presence_rising.tsv``
    -> ``("5mo_16mo-40_AL", "single_sample_tTest")``.  Anything else raises.
    """
    m = _NAME.match(os.path.basename(path))
    if not m:
        raise ValueError(f"not a baseline-presence rising table: {path}")
    return m.group("comparison"), m.group("test")


def normalise_origin(origin: str, earlier: str, later: str) -> str:
    """The rising table's ``origin_any_mouse`` without timepoint names.

    Example (5mo -> 16mo): ``de_novo_candidate_below_detection_at_5mo`` ->
    ``de_novo_below_detection``; ``5mo_not_covered`` -> ``baseline_not_covered``.
    An unknown label raises (a label from another comparison means the file
    name and its contents disagree).
    """
    mapping = {
        "standing_variation": "standing_variation",
        "de_novo_candidate": "de_novo",
        f"de_novo_candidate_below_detection_at_{earlier}": "de_novo_below_detection",
        f"allele_not_present_at_{later}": "not_present_later",
        f"{earlier}_not_covered": "baseline_not_covered",
    }
    if origin not in mapping:
        raise ValueError(
            f"unknown origin_any_mouse {origin!r} for {earlier} -> {later}"
        )
    return mapping[origin]


def load_rising(rising_dir: str) -> pd.DataFrame:
    """Every rising table in ``rising_dir``, with comparison / test / ages added.

    Header-only tables (no significant sites, or none that rose) add no rows
    and are logged.  No rising table at all raises; so does a missing column.
    Group names are read as text ("40" stays "40").
    """
    paths = sorted(
        glob.glob(os.path.join(rising_dir, "*_baseline_presence_rising.tsv"))
    )
    if not paths:
        raise FileNotFoundError(f"no *_baseline_presence_rising.tsv in {rising_dir}")
    text = {
        c: str
        for c in ("mag_id", "contig", "gene_id", "group_analyzed", "allele", "rose_in")
    }
    frames = []
    for path in paths:
        comparison, test = parse_rising_name(path)
        earlier, later, _, _ = parse_comparison(comparison)
        df = pd.read_csv(path, sep="\t", dtype=text)
        missing = [c for c in REQUIRED if c not in df.columns]
        if missing:
            raise ValueError(f"{os.path.basename(path)} lacks columns {missing}")
        if df.empty:
            logger.info(f"{comparison} {test}: no rising alleles (header-only table)")
            continue
        # The carriers (mice with the allele present at the later timepoint),
        # written by baseline presence under the comparison's own later label;
        # required like every other column (a table without it is from an older
        # baseline presence and must be regenerated).
        carriers_col = f"{later}_mice_allele_present"
        if carriers_col not in df.columns:
            raise ValueError(
                f"{os.path.basename(path)} lacks column {carriers_col} (the carriers): "
                f"rerun alleleflux-baseline-presence"
            )
        carriers = df[carriers_col].fillna("").astype(str)
        df = df[REQUIRED].copy()
        df["carriers"] = carriers
        df["group_analyzed"] = df["group_analyzed"].fillna("")  # blank = between-group
        df["gene_id"] = df["gene_id"].fillna("")  # intergenic
        df.insert(0, "age", later)
        df.insert(0, "period", f"{earlier}_{later}")
        df.insert(0, "test_type", test)
        df.insert(0, "comparison", comparison)
        df["origin"] = [
            normalise_origin(o, earlier, later) for o in df["origin_any_mouse"]
        ]
        frames.append(df)
    if not frames:
        return pd.DataFrame(
            columns=["comparison", "test_type", "period", "age"]
            + REQUIRED
            + ["carriers", "origin"]
        )
    return pd.concat(frames, ignore_index=True)


def load_classification(
    table: pd.DataFrame, metric: str, status_column: str
) -> tuple[int, pd.DataFrame]:
    """One floor's verdicts: ``(floor, [mag_id, group, period, status])``.

    ``table`` is a raw replacement_classification table.  Keeps the rows of
    ``metric`` (``dominant_strain_change`` = conANI), checks the contract --
    one floor (``min_voters``) per table, one verdict per MAG x group x period,
    only known statuses -- and renames to this script's keys.

    Example: a floor-8 table's row (MRGM_0016, 1D, 5mo_16mo, replaced) ->
    ``(8, [MRGM_0016, "1D", "5mo_16mo", "replaced"])``.
    """
    t = table[table["metric"] == metric].copy()
    if t.empty:
        raise ValueError(f"no rows for metric {metric!r} in the classification table")
    floors = t["min_voters"].unique()
    if len(floors) != 1:
        raise ValueError(
            f"one table must hold one voter floor, found min_voters = {sorted(floors)}"
        )
    key = ["MAG_ID", "group", "transition"]
    dup = t.duplicated(key)
    if dup.any():
        raise ValueError(
            f"{int(dup.sum())} duplicate MAG x group x period rows, e.g. "
            f"{t.loc[dup, key].iloc[0].tolist()}"
        )
    unknown = set(t[status_column]) - set(STATUS_TO_CLASS)
    if unknown:
        raise ValueError(
            f"unknown strain statuses in {status_column}: {sorted(unknown)}"
        )
    t = t.rename(
        columns={
            "MAG_ID": "mag_id",
            "transition": "period",
            status_column: "status",
        }
    )
    t["group"] = t["group"].astype(str)  # "40" must never become the integer 40
    t["mag_id"] = t["mag_id"].astype(str)
    return int(floors[0]), t[["mag_id", "group", "period", "status"]]


def combine_arms(status_a: str, status_b: str) -> str:
    """Class of a between-group allele from its comparison's two arms.

    Rules, in order: any arm replaced -> ``replaced``; both kept ->
    ``not_replaced``; any arm tied -> ``tie``; otherwise ``unverified``.

    >>> combine_arms("not_replaced", "replaced")
    'replaced'
    >>> combine_arms("not_replaced", "no_voters")
    'unverified'
    """
    classes = []
    for status in (status_a, status_b):
        if status not in STATUS_TO_CLASS:
            raise ValueError(
                f"unknown strain status {status!r}; expected one of {sorted(STATUS_TO_CLASS)}"
            )
        classes.append(STATUS_TO_CLASS[status])
    if "replaced" in classes:
        return "replaced"  # one swapped arm is enough to explain the difference
    if classes == ["not_replaced", "not_replaced"]:
        return "not_replaced"
    if "tie" in classes:
        return "tie"
    return "unverified"


def _status_of(
    rows: pd.DataFrame, verdicts: pd.DataFrame, group_col: str
) -> np.ndarray:
    """Each row's raw verdict for MAG x ``group_col`` x period; missing -> not_in_strain_table."""
    v = verdicts.rename(columns={"group": group_col})
    merged = rows[["mag_id", group_col, "period"]].merge(
        v, on=["mag_id", group_col, "period"], how="left", validate="many_to_one"
    )
    return merged["status"].fillna("not_in_strain_table").to_numpy()


def attach_strain(
    rising: pd.DataFrame, tables: list[tuple[int, pd.DataFrame]]
) -> pd.DataFrame:
    """One copy of ``rising`` per floor, with ``strain_status`` and ``strain_class``.

    ``tables`` = ``load_classification`` results, one per floor.  Within-group
    rows read the arm they rose in; between-group rows read both arms of their
    comparison and combine them (``strain_status`` then shows both, e.g.
    ``40:not_replaced|AL:replaced``).  REPORTING ONLY: no origin changes here;
    ``apply_strain_gate`` decides what the verdict does to a de novo allele.

    Example (floor 8, MAG_B replaced in 40): MAG_B's C gets strain_status
    ``replaced``, class ``replaced``; at floor 4, where MAG_B is kept,
    ``not_replaced``.  Its origin is ``de_novo_below_detection`` in both copies.
    """
    floors = [f for f, _ in tables]
    if len(set(floors)) != len(floors):
        raise ValueError(f"two classification tables share a voter floor: {floors}")
    copies = []
    for floor, verdicts in tables:
        out = rising.copy()
        out.insert(0, "floor", floor)
        is_within = (out["group_analyzed"] != "").to_numpy()
        status = np.empty(len(out), dtype=object)
        cls = np.empty(len(out), dtype=object)
        if is_within.any():
            raw = _status_of(out[is_within], verdicts, "rose_in")
            status[is_within] = raw
            cls[is_within] = [STATUS_TO_CLASS[s] for s in raw]
        if (~is_within).any():
            between = out[~is_within].copy()
            groups = [parse_comparison(c)[2:] for c in between["comparison"]]
            between["arm_a"] = [g[0] for g in groups]
            between["arm_b"] = [g[1] for g in groups]
            a = _status_of(between, verdicts, "arm_a")
            b = _status_of(between, verdicts, "arm_b")
            status[~is_within] = [
                f"{ga}:{sa}|{gb}:{sb}"
                for ga, sa, gb, sb in zip(between["arm_a"], a, between["arm_b"], b)
            ]
            cls[~is_within] = [combine_arms(sa, sb) for sa, sb in zip(a, b)]
        out["strain_status"] = status
        out["strain_class"] = cls
        copies.append(out)
    return pd.concat(copies, ignore_index=True)


def attach_carrier_gate(
    alleles: pd.DataFrame, turnover_dir: str, min_determined: int
) -> pd.DataFrame:
    """Add the carrier-level strain verdict to every row (spec 2.4).

    Reads ``{turnover_dir}/{mag}_strain_turnover.tsv`` once per MAG (the
    per-mouse verdicts of alleleflux-strain-turnover) and asks
    ``carrier_strain_verdicts`` about each row's carriers for its period.
    Adds ``carrier_strain_status``, ``n_carriers``, ``n_carriers_determined``
    and ``n_carriers_strain_changed``; changes no origin.  A MAG with no
    turnover table gives ``strain_undetermined`` for all its rows, with one
    warning per MAG.

    Example: carriers m1, m2, m3 with verdicts True, True, False at floor 3 ->
    ``new_strain_carriers``, 3, 3, 2.
    """
    out = alleles.copy()
    tables: dict[str, pd.DataFrame | None] = {}
    for mag in sorted(out["mag_id"].unique()):
        path = os.path.join(turnover_dir, f"{mag}_strain_turnover.tsv")
        if not os.path.exists(path):
            logger.warning(
                f"{mag}: no strain-turnover table at {path}; its alleles are "
                f"strain_undetermined under the carrier gate"
            )
            tables[mag] = None
            continue
        tables[mag] = pd.read_csv(
            path,
            sep="\t",
            dtype={"subjectID": str, "transition": str, "sample_t2": str},
            usecols=["subjectID", "transition", "sample_t2", "dominant_strain_change"],
        )
    status, n_car, n_det, n_chg = [], [], [], []
    for mag, period, carriers in zip(out["mag_id"], out["period"], out["carriers"]):
        mice = [m for m in str(carriers).split(",") if m]
        table = tables[mag]
        if table is None:
            verdict = (CARRIER_UNDETERMINED, 0, 0)
        else:
            verdict = carrier_strain_verdicts(mice, table, period, min_determined)
        status.append(verdict[0])
        n_car.append(len(mice))
        n_det.append(verdict[1])
        n_chg.append(verdict[2])
    out["carrier_strain_status"] = status
    out["n_carriers"] = n_car
    out["n_carriers_determined"] = n_det
    out["n_carriers_strain_changed"] = n_chg
    logger.info(
        f"carrier gate (floor {min_determined} determined carriers): "
        f"{pd.Series(status).value_counts().to_dict()}"
    )
    return out


def _group_verdict(strain_class: str):
    """What the group vote does to a de novo allele (user ruling 2026-10-23):
    replaced -> new_strain; not_replaced -> None (kept); anything else (unverified:
    too few / no voters, not in the table; or tie) -> strain_undetermined."""
    if strain_class == "replaced":
        return "new_strain"
    if strain_class == "not_replaced":
        return None
    return CARRIER_UNDETERMINED


def _carrier_verdict(carrier_status: str):
    """What the carriers' verdict does to a de novo allele."""
    if carrier_status == CARRIER_NEW_STRAIN:
        return "new_strain"
    if carrier_status == CARRIER_SAME_STRAIN:
        return None
    if carrier_status == CARRIER_UNDETERMINED:
        return CARRIER_UNDETERMINED
    raise ValueError(f"unknown carrier_strain_status {carrier_status!r}")


def apply_strain_gate(alleles: pd.DataFrame, gate: str) -> pd.DataFrame:
    """Which strain verdict decides whether a de novo allele keeps its origin (spec 2.4).

    The ONLY place an origin changes after the cluster rule, and only for the
    two de novo flavours (``DE_NOVO``); standing alleles, block events and the
    non-verdicts are never touched.  ``origin_any_mouse`` keeps the raw verdict.

    BOTH verdict columns stay on every row whatever the gate: ``strain_class``
    (group) and ``carrier_strain_status`` (carrier) are written whenever their
    input was given, so the two can be compared on the same run.  The gate only
    decides which of them is READ; under ``group`` the carrier columns are
    reporting only, under ``carrier`` the group columns are.

    * ``group``: from ``strain_class`` (``attach_strain``): replaced ->
      ``new_strain``; not_replaced -> kept; anything else ->
      ``strain_undetermined``.
    * ``carrier``: from ``carrier_strain_status`` (``attach_carrier_gate``):
      new_strain_carriers -> ``new_strain``; same_strain_carriers -> kept;
      strain_undetermined -> ``strain_undetermined``.
    * ``both``: an allele must pass both.  If either gate is undetermined the
      allele is ``strain_undetermined`` (removed, not explained: one gate
      could not speak); else if either says a swap, ``new_strain``; else kept
      (``strain_undetermined`` over ``new_strain`` over kept).
    * ``none``: no relabel.

    A gate whose column is missing raises (the inputs were not given).

    Example: group not_replaced, carriers new_strain_carriers -> ``group`` keeps
    de_novo, ``carrier`` and ``both`` give new_strain.  Group replaced,
    carriers undetermined -> ``group`` gives new_strain, ``carrier`` and
    ``both`` give strain_undetermined.
    """
    if gate not in STRAIN_GATES:
        raise ValueError(
            f"unknown strain gate {gate!r}; expected one of {STRAIN_GATES}"
        )
    out = alleles.copy()
    if gate == "none":
        return out
    need = {"group": ["strain_class"], "carrier": ["carrier_strain_status"]}
    need["both"] = need["group"] + need["carrier"]
    missing = [c for c in need[gate] if c not in out.columns]
    if missing:
        raise ValueError(
            f"--strain_gate {gate} needs {missing} on every row: give "
            f"--classification (group) / --turnover_dir (carrier)"
        )
    is_dn = out["origin"].isin(DE_NOVO).to_numpy()
    verdicts = []
    for i in np.flatnonzero(is_dn):
        row = out.iloc[i]
        candidates = []
        if gate in ("group", "both"):
            candidates.append(_group_verdict(row["strain_class"]))
        if gate in ("carrier", "both"):
            candidates.append(_carrier_verdict(row["carrier_strain_status"]))
        verdicts.append(max(candidates, key=_GATE_RANK.__getitem__))
    idx = out.index[is_dn]
    new = pd.Series(verdicts, index=idx, dtype=object)
    out.loc[idx[new.notna().to_numpy()], "origin"] = new.dropna()
    return out


def summarise(alleles: pd.DataFrame) -> pd.DataFrame:
    """Counts and percentages per floor x comparison x group x strain class x origin.

    Without strain columns: ``strain_class = "all"`` and ``floor`` NaN.  With
    them: each strain class, plus ``all`` (every class of that floor x
    comparison x group together; floors and comparisons are never pooled).
    Every origin gets a row, zero counts included, so a plot never has gaps.
    ``pct_of_all`` divides by ``n_rising`` (all rising alleles in the block),
    ``pct_of_called`` by ``n_called`` (the four verdicts) and is NaN for the two
    non-verdicts.
    """
    gated = (
        "strain_class" in alleles.columns or "carrier_strain_status" in alleles.columns
    )
    if "strain_class" in alleles.columns:
        tagged = pd.concat(
            [alleles, alleles.assign(strain_class="all")], ignore_index=True
        )
    else:
        tagged = alleles.assign(strain_class="all", floor=np.nan)
    # Without any strain input the two gate outcomes cannot occur: no zero rows for them.
    origins = (
        ORIGINS
        if gated
        else [o for o in ORIGINS if o not in ("new_strain", CARRIER_UNDETERMINED)]
    )
    rows = []
    for block, g in tagged.groupby(BLOCK, sort=True, dropna=False):
        n_rising = len(g)
        n_called = int(g["origin"].isin(CALLED).sum())
        for origin in origins:
            hit = g[g["origin"] == origin]
            n = len(hit)
            rows.append(
                {
                    **dict(zip(BLOCK, block)),
                    "origin": origin,
                    "n_alleles": n,
                    "n_mags": hit["mag_id"].nunique(),
                    "n_rising": n_rising,
                    "n_called": n_called,
                    "pct_of_all": 100 * n / n_rising if n_rising else np.nan,
                    "pct_of_called": (
                        100 * n / n_called if origin in CALLED and n_called else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def tally_the_origins(args: argparse.Namespace) -> int:
    """Orchestrator: rising tables (+ classification tables) -> allele table + summary."""
    rising = load_rising(args.rising_dir)
    logger.info(
        f"{len(rising):,} rising alleles from {rising['comparison'].nunique()} comparisons"
    )
    # Cluster rule first (spec 2.3), on every candidate, before any strain gate.
    rising = assign_blocks(rising, args.cluster_window_bp, args.max_de_novo_neighbours)
    n_block_rows = int((rising["origin"] == BLOCK_EVENT).sum())
    n_blocks = rising.loc[rising["block_id"] != "", "block_id"].nunique()
    logger.info(
        f"cluster rule (window {args.cluster_window_bp} bp, > "
        f"{args.max_de_novo_neighbours} de novo neighbours): {n_block_rows:,} de novo "
        f"rows relabelled block_event in {n_blocks} blocks"
    )
    if args.classification:
        tables = [
            load_classification(
                pd.read_csv(path, sep="\t", dtype={"group": str, "MAG_ID": str}),
                args.strain_metric,
                args.status_column,
            )
            for path in args.classification
        ]
        alleles = attach_strain(rising, tables)
        logger.info(
            f"group strain verdicts attached: floors {sorted(f for f, _ in tables)}"
        )
    else:
        alleles = rising
    if args.turnover_dir:
        alleles = attach_carrier_gate(
            alleles, args.turnover_dir, args.min_determined_carriers
        )
    # The gate: the only place an origin changes after the cluster rule.
    alleles = apply_strain_gate(alleles, args.strain_gate)
    logger.info(
        f"strain gate {args.strain_gate}: "
        f"{int((alleles['origin'] == 'new_strain').sum()):,} de novo rows -> new_strain, "
        f"{int((alleles['origin'] == CARRIER_UNDETERMINED).sum()):,} -> strain_undetermined"
    )
    summary = summarise(alleles)
    os.makedirs(args.output_dir, exist_ok=True)
    allele_path = os.path.join(args.output_dir, "rising_origin_alleles.tsv")
    summary_path = os.path.join(args.output_dir, "rising_origin_summary.tsv")
    alleles.to_csv(allele_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)
    logger.info(f"wrote {len(alleles):,} allele rows to {allele_path}")
    logger.info(f"wrote {len(summary):,} summary rows to {summary_path}")
    return 0


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Standing variation or de novo for the alleles that rose: "
        "one tidy summary over every comparison's rising table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--rising_dir",
        required=True,
        help="folder of *_baseline_presence_rising.tsv tables (baseline_presence/)",
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--classification",
        nargs="*",
        default=None,
        help="replacement_classification.tsv table(s), one per voter floor; when "
        "given, de novo in a replaced MAG becomes new_strain and the summary is "
        "split by floor and strain class",
    )
    parser.add_argument(
        "--cluster_window_bp",
        type=int,
        default=1000,
        help="cluster rule: two de novo alleles on the same contig that rose in the "
        "same group are neighbours when at most this many bp apart (any age); "
        "neighbours chain into blocks",
    )
    parser.add_argument(
        "--max_de_novo_neighbours",
        type=int,
        default=0,
        help="cluster rule: a de novo allele with MORE than this many de novo "
        "neighbours becomes a block_event (one imported or recombined chunk), never "
        "an independent de novo mutation; 0 = only lone alleles stay de novo",
    )
    parser.add_argument(
        "--turnover_dir",
        default=None,
        help="folder of {mag}_strain_turnover.tsv tables (alleleflux-strain-turnover): "
        "per-mouse strain verdicts the carrier gate reads for each allele's carriers",
    )
    parser.add_argument(
        "--strain_gate",
        choices=STRAIN_GATES,
        default="carrier",
        help="which strain verdict decides whether a de novo allele stands: carrier = "
        "the allele's own carriers (needs --turnover_dir); group = the across-cage "
        "vote (needs --classification); both = the stricter of the two; none = no "
        "strain gate.  Both verdicts are written to every row whenever given",
    )
    parser.add_argument(
        "--min_determined_carriers",
        type=int,
        default=3,
        help="carrier gate: carriers with a per-mouse verdict needed before the "
        "carriers' strain is called; fewer -> strain_undetermined",
    )
    parser.add_argument(
        "--strain_metric",
        default="dominant_strain_change",
        help="classification metric: dominant_strain_change (conANI) or "
        "strain_replacement (popANI, stricter)",
    )
    parser.add_argument(
        "--status_column",
        default="strain_status_replicates",
        help="verdict column: strain_status_replicates (replicates vote) or "
        "strain_status (mice vote)",
    )
    args = parser.parse_args()
    if args.strain_gate in ("carrier", "both") and not args.turnover_dir:
        parser.error(f"--strain_gate {args.strain_gate} needs --turnover_dir")
    if args.strain_gate in ("group", "both") and not args.classification:
        parser.error(f"--strain_gate {args.strain_gate} needs --classification")
    return tally_the_origins(args)


if __name__ == "__main__":
    raise SystemExit(main())
