"""Gene-level de novo accumulation: which genes collect new alleles, per cage.

The question (spec section 3).  Every cage keeps a diary of the new letters
that appeared in a genome between the two timepoints of a comparison.  Each
entry is filed under the gene it landed in.  Test A asks whether a gene got
more entries than its share of the readable genome inside one group; Test B
asks whether a gene collected entries per readable base faster in one group
than in the other.  The cage is the replicate throughout.

This module holds the COUNTING side as pure functions over frames (the CLI
streams the real files through them; the two Poisson tests live in
``statistics/poisson_rate_tests.py``):

* ``cage_rises``: from the pipeline's mean-change table (``allele_freq.py``,
  one row per position x cage: the cage-average change in each base's share,
  over the mice sampled at both timepoints, after the zero-diff filter) the
  letters whose share went UP, one row per (cage, position, letter).  These
  are the candidates; nothing else can be a de novo event;
* ``paired_samples``: the samples of the mice sampled at both timepoints --
  the ONLY mice the mean-change table saw, hence the only reads that count
  as "searched" for the denominator (a cage of unpaired mice is "not
  searched", never "no mutations");
* ``CoverageMask``: per cage and timepoint, one bit per contig position, set
  where some paired sample has at least ``min_cov`` reads (the same
  ``min_cov`` the presence verdicts use).  Readable = set at both timepoints.
  A cage is the pair (group, replicate): replicate labels may repeat across
  groups (the bundled data names the cages of both groups A..D), so the
  replicate alone would merge two cages;
  Bits, not rows: a genome of 5 Mb and 100 cages is 1 GB of booleans, where
  a frame of readable positions would be 500 million rows;
* ``covered_bp_by_gene`` / ``covered_bp_total``: readable positions inside a
  gene's span (Prodigal coordinates made 0-based inclusive, as
  ``profile_mags.map_genes`` does) and over the whole genome -- the
  denominators of the two tests.  Readable bases, never gene length: a gene
  read over 300 of its 1,200 bases was searched over 300.

Worked cage (min_cov 5).  Cage A has mice m1 and m2 at both timepoints.
Position 10: m1 has 5 reads before and 4 after, m2 has 7 after -> readable
(some paired mouse >= 5 at each timepoint).  Position 11: m1 has 9 reads
before and no row after -> not readable.  Gene g1 spans 100..200 and cage A
is readable at 99, 100, 150, 201 -> ``covered_bp`` 2, genome total 4.  The
mean-change table says at position 10 the share of T rose by 0.30 in cage A
-> one candidate row (cageA, c1, 10, T).

Positions are 0-based everywhere (cache and ``map_genes``); a position at or
past the contig length is an error, never silently dropped.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from Bio.Seq import Seq

from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
    NUCLEOTIDES,
)
from alleleflux.scripts.analysis.baseline_presence import (
    EVIDENCE_BELOW_DETECTION,
    EVIDENCE_NOT_COVERED,
    Timepoints,
    assess_rows,
    baseline_cutoff,
    baseline_rows,
)
from alleleflux.scripts.analysis.de_novo_rules import (
    BLOCK_EVENT,
    CARRIER_NEW_STRAIN,
    CARRIER_UNDETERMINED,
    DE_NOVO_ORIGINS,
    assign_blocks,
    carrier_strain_verdicts,
)

# load_genes is re-exported here: tests and callers import it from this module.
from alleleflux.scripts.utilities.genes import (  # noqa: F401 (re-export)
    COMPLEMENT_MAP,
    coding_index,
    load_genes,
)

logger = logging.getLogger(__name__)

# The per-cage origin words.  The three "counted" ones come from
# de_novo_rules (DE_NOVO_ORIGINS + BLOCK_EVENT); the rest mirror rising
# origin's vocabulary so a reader of both tables meets one set of words.
ORIGIN_STANDING = "standing_variation"
ORIGIN_DE_NOVO, ORIGIN_DE_NOVO_BELOW_DETECTION = DE_NOVO_ORIGINS
ORIGIN_NEW_STRAIN = "new_strain"
ORIGIN_NOT_PRESENT_LATER = "not_present_later"
ORIGIN_BASELINE_NOT_COVERED = "baseline_not_covered"
BASELINE_SCOPES = ("individual", "own_group", "comparison", "cohort")
STRAIN_GATES = ("carrier", "none")  # no group vote at cage level
BLOCK_RULES = ("exclude", "once")

VERDICT_COLUMNS = [
    "origin",
    "block_id",
    "n_de_novo_neighbours",
    "carriers",
    "n_carriers",
    "carrier_strain_status",
    "n_carriers_determined",
    "n_carriers_strain_changed",
    "n_baseline_samples_covered",
    "n_baseline_samples_present",
]

# The mean-change table's per-base columns (``A_frequency_diff_mean`` ...) and
# the one-letter allele each stands for.  Derived from the ONE constant the
# allele-frequency step uses for its frequency columns, so a renamed column
# there breaks loudly here instead of silently matching nothing.
DIFF_MEAN_COLUMNS = {f"{col}_diff_mean": col.split("_")[0] for col in NUCLEOTIDES}

RISE_COLUMNS = [
    "contig",
    "position",
    "gene_id",
    "replicate",
    "group",
    "allele",
    "diff_mean",
    "subjectID_count",
]


def cage_rises(mean_changes: pd.DataFrame, groups: tuple[str, str]) -> pd.DataFrame:
    """The letters whose cage-average share went up, one row per (cage, position, letter).

    ``mean_changes`` is one MAG's ``{mag}_allele_frequency_changes_mean.parquet``
    (columns ``contig, gene_id, position, replicate, group,
    {A,T,G,C}_frequency_diff_mean, subjectID_count``; the diff is later minus
    earlier in the comparison's listed order, averaged over the cage's paired
    mice).  ``groups`` are the comparison's two group names; rows of any
    other group are dropped.  Labels come back as plain text (the parquet
    stores categoricals whose category order differs between runs) and a
    missing ``gene_id`` becomes ``""`` so "between genes" has one spelling.

    A rise is strictly ``> 0``: a letter whose share did not move is not a
    candidate, and a letter that fell is the other side of someone else's
    rise.  At a three-letter site two letters can rise at once; both become
    rows (the de novo checks judge each letter on its own).

    Example: cage A at position 10 with A -0.3 / T +0.3 and cage B at 30 with
    A +0.1 / T -0.1 -> two rows: (cageA, 10, T, 0.3) and (cageB, 30, A, 0.1).
    """
    missing = [c for c in DIFF_MEAN_COLUMNS if c not in mean_changes.columns]
    if missing:
        raise ValueError(
            f"mean-change table lacks the per-base change column(s) {missing}; "
            f"expected {sorted(DIFF_MEAN_COLUMNS)}"
        )
    wanted = {str(g) for g in groups}
    df = mean_changes.copy()
    # Text, not categories: a category code is meaningless across tables.
    df["replicate"] = df["replicate"].astype(str)
    df["group"] = df["group"].astype(str)
    df = df[df["group"].isin(wanted)]
    df["gene_id"] = df["gene_id"].fillna("").astype(str)
    df["contig"] = df["contig"].astype(str)
    df["position"] = df["position"].astype("int64")
    # One stacked frame of (row, letter, change), then keep the rises.
    pieces = []
    for col, allele in DIFF_MEAN_COLUMNS.items():
        change = df[col].astype("float64")
        up = change > 0
        if not up.any():
            continue
        piece = df.loc[
            up, ["contig", "position", "gene_id", "replicate", "group"]
        ].copy()
        piece["allele"] = allele
        piece["diff_mean"] = change[up].to_numpy()
        piece["subjectID_count"] = (
            df.loc[up, "subjectID_count"].astype("int64").to_numpy()
        )
        pieces.append(piece)
    if not pieces:
        return pd.DataFrame(columns=RISE_COLUMNS)
    out = pd.concat(pieces, ignore_index=True)
    # Deterministic order: by cage, contig, position, then letter.
    out = out.sort_values(["replicate", "contig", "position", "allele"], kind="stable")
    return out[RISE_COLUMNS].reset_index(drop=True)


def paired_samples(sheet: pd.DataFrame, earlier: str, later: str) -> pd.DataFrame:
    """The samples of mice present at BOTH timepoints: the denominator's mice.

    The mean-change table is built from exactly these mice (``allele_freq.py``
    inner-joins each subject's two timepoints), so only their reads can have
    revealed a rise; a cage made of unpaired mice was never searched, and its
    readable bases must be zero rather than "zero mutations in N bases".
    ``sheet`` is the run's metadata (``sample_id, subjectID, replicate, group,
    time``; other columns ignored); labels are compared as text.  A mouse with
    two samples at one of the timepoints raises, naming it: these designs
    sample a mouse once per timepoint, and picking one would hide a sheet
    error.

    Example: m1 (pre, post), m2 (pre only), m4 (pre, mid, post) -> the four
    samples of m1 and m4 at pre and post; m2's and the mid sample are out.
    """
    cols = ["sample_id", "subjectID", "replicate", "group", "time"]
    rows = sheet[cols].astype(str)
    rows = rows[rows["time"].isin([str(earlier), str(later)])]
    counts = rows.groupby(["subjectID", "time"]).size()
    doubled = counts[counts > 1]
    if len(doubled):
        who = sorted({s for s, _ in doubled.index})
        raise ValueError(
            f"mouse/mice {who} have more than one sample at a timepoint of "
            f"[{earlier}, {later}]; one sample per mouse and timepoint is required"
        )
    per_subject = rows.groupby("subjectID")["time"].nunique()
    both = per_subject[per_subject == 2].index
    return rows[rows["subjectID"].isin(both)].reset_index(drop=True)


class CoverageMask:
    """One bit per contig position, per cage and timepoint: "some sample had >= min_cov reads here".

    Built up from cache batches (``cache_io.iter_cache_batches``; the caller
    restricts the samples to ``paired_samples``), then asked for the readable
    positions of a comparison: set at ``earlier`` AND at ``later``.

    ``contig_lengths`` (``utilities.build_contig_length_index``) sizes the
    arrays; a row at a position >= its contig's length, or on an unknown
    contig, raises -- a silent drop would shrink the denominator.

    Example (min_cov 5): rows m1_pre c1:10 cov 5, m1_post c1:10 cov 4, m2_post
    c1:10 cov 7, m1_pre c1:11 cov 9 (cage A) -> readable("pre", "post")
    ["cageA"]["c1"] is True at 10 only (11 has no later-timepoint row).
    """

    def __init__(self, contig_lengths: dict[str, int], min_cov: int):
        self.contig_lengths = {str(c): int(n) for c, n in contig_lengths.items()}
        self.min_cov = int(min_cov)
        # (group, replicate, timepoint) -> {contig -> bool array}; arrays are
        # created on first sight of a (cage, timepoint) so idle cages cost nothing.
        self._bits: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}

    def _arrays(
        self, group: str, replicate: str, timepoint: str
    ) -> dict[str, np.ndarray]:
        key = (group, replicate, timepoint)
        if key not in self._bits:
            self._bits[key] = {
                c: np.zeros(n, dtype=bool) for c, n in self.contig_lengths.items()
            }
        return self._bits[key]

    def add(self, rows: pd.DataFrame) -> None:
        """Set the bits for one batch of cache rows.

        ``rows`` need ``group, replicate, timepoint, contig, position,
        total_coverage`` (``sample_id`` is not consulted: which samples are in
        is the caller's filter).  Rows under ``min_cov`` are ignored; the rest
        set ``bits[(group, replicate, timepoint)][contig][position]``.
        """
        deep = rows[rows["total_coverage"].astype("int64") >= self.min_cov]
        if deep.empty:
            return
        unknown = sorted(set(deep["contig"].astype(str)) - set(self.contig_lengths))
        if unknown:
            raise ValueError(
                f"cache rows on contig(s) {unknown} that the contig length index "
                "does not know; the MAG mapping and the cache disagree"
            )
        for (group, replicate, timepoint, contig), part in deep.groupby(
            [
                deep["group"].astype(str),
                deep["replicate"].astype(str),
                deep["timepoint"].astype(str),
                deep["contig"].astype(str),
            ],
            sort=False,
        ):
            positions = part["position"].to_numpy(dtype="int64")
            length = self.contig_lengths[contig]
            if positions.size and (positions.max() >= length or positions.min() < 0):
                raise ValueError(
                    f"contig {contig}: position {int(positions.max())} outside "
                    f"0..{length - 1} (positions are 0-based; is the length index stale?)"
                )
            self._arrays(group, replicate, timepoint)[contig][positions] = True

    def cages(self) -> list[tuple[str, str]]:
        """Every (group, replicate) seen at any timepoint, sorted."""
        return sorted({(group, replicate) for group, replicate, _ in self._bits})

    def readable(
        self, earlier: str, later: str
    ) -> dict[tuple[str, str], dict[str, np.ndarray]]:
        """Per cage (group, replicate), per contig: True where deep at BOTH timepoints.

        A cage seen at only one timepoint gets all-False arrays (present in
        the result so it appears in the tables with zeros, never missing).
        """
        out: dict[tuple[str, str], dict[str, np.ndarray]] = {}
        for group, replicate in self.cages():
            t0 = self._bits.get((group, replicate, str(earlier)))
            t1 = self._bits.get((group, replicate, str(later)))
            if t0 is None or t1 is None:
                out[(group, replicate)] = {
                    c: np.zeros(n, dtype=bool) for c, n in self.contig_lengths.items()
                }
                continue
            out[(group, replicate)] = {c: t0[c] & t1[c] for c in self.contig_lengths}
        return out


def covered_bp_by_gene(
    readable: dict[tuple[str, str], dict[str, np.ndarray]], genes: pd.DataFrame
) -> pd.DataFrame:
    """Readable positions inside each gene, per cage: one row per (cage, gene), zeros kept.

    ``genes`` has ``gene_id, contig, start, end`` with 0-based INCLUSIVE
    coordinates (``profile_mags.map_genes`` subtracts 1 from both Prodigal
    ends).  A gene on a contig the mask does not know raises: the gene file
    and the reference disagree.

    Example: cage A readable at 99, 100, 150, 201 on c1; gene g1 = 100..200
    -> ``covered_bp`` 2 (100 and 150; 99 and 201 are outside).
    """
    records = []
    gene_rows = list(
        zip(
            genes["gene_id"].astype(str),
            genes["contig"].astype(str),
            genes["start"].astype("int64"),
            genes["end"].astype("int64"),
        )
    )
    for group, replicate in sorted(readable):
        arrays = readable[(group, replicate)]
        for gene_id, contig, start, end in gene_rows:
            if contig not in arrays:
                raise ValueError(
                    f"gene {gene_id} sits on contig {contig}, which the coverage "
                    "mask does not know; the gene file and the reference disagree"
                )
            # end is inclusive, hence the + 1; a span past the contig end
            # simply counts what exists (numpy slices clip).
            covered = int(arrays[contig][start : end + 1].sum())
            records.append(
                {
                    "group": group,
                    "replicate": replicate,
                    "gene_id": gene_id,
                    "covered_bp": covered,
                }
            )
    return pd.DataFrame(
        records, columns=["group", "replicate", "gene_id", "covered_bp"]
    )


def covered_bp_total(
    readable: dict[tuple[str, str], dict[str, np.ndarray]],
) -> pd.DataFrame:
    """Readable positions over the whole genome, per cage (the genome-wide denominator)."""
    records = [
        {
            "group": group,
            "replicate": replicate,
            "covered_bp_total": int(sum(int(a.sum()) for a in arrays.values())),
        }
        for (group, replicate), arrays in sorted(readable.items())
    ]
    return pd.DataFrame(records, columns=["group", "replicate", "covered_bp_total"])


# ---------------------------------------------------------------------------
# Task 2: the Phase 2 rules, cage by cage
# ---------------------------------------------------------------------------


def _cage_key(group: str, replicate: str) -> str:
    """The cluster key handed to ``assign_blocks`` as ``rose_in``: one cage."""
    return f"{group}|{replicate}"


def classify_cage_rises(
    rises: pd.DataFrame,
    cache_rows: pd.DataFrame,
    sheet: pd.DataFrame,
    *,
    mag_id: str,
    tps: Timepoints,
    scope: str,
    model: np.ndarray,
    min_freq: float,
    min_cov: int,
    window_bp: int,
    max_neighbours: int,
    strain_gate: str,
    min_determined_carriers: int,
    turnover: pd.DataFrame | None,
) -> pd.DataFrame:
    """Phase 2's three questions, asked of every candidate rise with the CAGE as the unit.

    For each row of ``rises`` (one MAG; Task 1's ``cage_rises``) the cache rows
    of that site for EVERY loaded sample (``cache_rows``: ``contig, position,
    sample_id, subjectID, replicate, group, timepoint, A, C, G, T``; under the
    cohort scope every timepoint of the run, otherwise the comparison's two)
    are judged in this fixed order, reusing the Phase 2 functions as they are:

    1. **presence verdicts** (``baseline_presence.assess_rows``): per sample,
       ``present`` / ``below_detection`` / ``absent`` / ``not_covered`` for
       the rising letter, with the bar (``model``), ``min_freq`` and ``min_cov``;
    2. **carriers**: the cage's later-timepoint mice with the letter
       ``present``.  None -> ``not_present_later`` (a cage-average tick with no
       mouse clearing the bar is noise);
    3. **baseline** (``baseline_presence.baseline_rows`` with ``scope``):
       ``individual`` = the carriers' own earlier samples, ``own_group`` = the
       cage's group's earlier samples, ``comparison`` = both groups' earlier
       samples, ``cohort`` = every sample dated before the first carrier
       (``baseline_cutoff``).  No covered baseline sample ->
       ``baseline_not_covered``; any ``present`` -> ``standing_variation``; any
       ``below_detection`` -> ``de_novo_below_detection``; else ``de_novo``;
    4. **blocks** (``de_novo_rules.assign_blocks``) over the de novo rows of
       one cage (``rose_in`` = the cage): neighbours within ``window_bp`` on a
       contig, more than ``max_neighbours`` of them -> ``block_event`` with a
       ``block_id``.  Two cages never cluster together;
    5. **strain gate** on the two de novo flavours only:
       ``carrier`` asks ``de_novo_rules.carrier_strain_verdicts`` about the
       cage's carriers for the period ``{earlier}_{later}`` in ``turnover``
       (this MAG's ``{mag}_strain_turnover.tsv``; None = no table):
       ``new_strain_carriers`` -> ``new_strain``, undetermined ->
       ``strain_undetermined``, same strain -> kept.  ``none`` skips the gate.
       ``group`` / ``both`` are refused: a cage has no group vote.

    ``sheet`` supplies ``collection_date`` by ``sample_id``; without dates the
    cohort baseline is the comparison's earlier timepoint (``baseline_rows``).

    Returns ``rises``' columns plus ``VERDICT_COLUMNS``, same row order.
    Which rows count as events is ``is_event``'s decision, not this one's.

    Example (cohort, min_cov 5, bar 3 at 40x): cage A = m1 + m2, candidate T
    at 3000; every earlier sample has 40 A and 0 T, m1's later sample 20 A /
    20 T -> carriers "m1", 4 covered baseline samples, none present ->
    ``de_novo``; add 6 T in another group's earlier sample -> ``standing_variation``;
    add a second de novo rise in cage A at 3300 -> both ``block_event``.
    """
    if scope not in BASELINE_SCOPES:
        raise ValueError(
            f"baseline scope must be one of {BASELINE_SCOPES}, got {scope!r}"
        )
    if strain_gate not in STRAIN_GATES:
        raise ValueError(
            f"strain_gate must be one of {STRAIN_GATES} for the gene-level analysis "
            f"(a cage has no group vote), got {strain_gate!r}"
        )
    out_columns = list(rises.columns) + VERDICT_COLUMNS
    if rises.empty:
        return pd.DataFrame(columns=out_columns)

    # --- the rows every verdict is read from: labels as text, dates and the
    # sheet's ``time`` word attached by sample id (baseline_rows reads ``time``).
    rows = cache_rows.copy()
    for col in ("contig", "sample_id", "subjectID", "replicate", "group", "timepoint"):
        rows[col] = rows[col].astype(str)
    rows["position"] = rows["position"].astype("int64")
    rows["time"] = rows["timepoint"]
    meta = sheet.set_index(sheet["sample_id"].astype(str))
    rows["collection_date"] = (
        meta["collection_date"].reindex(rows["sample_id"]).to_numpy()
        if "collection_date" in meta.columns
        else np.nan
    )
    by_site = {
        key: idx for key, idx in rows.groupby(["contig", "position"]).indices.items()
    }

    # --- stages 1-3, one candidate at a time (the per-site frames are tiny:
    # one row per loaded sample).
    records = []
    for rise in rises.itertuples(index=False):
        contig, position = str(rise.contig), int(rise.position)
        group, replicate, allele = (
            str(rise.group),
            str(rise.replicate),
            str(rise.allele),
        )
        rec = {
            "carriers": "",
            "n_carriers": 0,
            "carrier_strain_status": "",
            "n_carriers_determined": 0,
            "n_carriers_strain_changed": 0,
            "n_baseline_samples_covered": 0,
            "n_baseline_samples_present": 0,
        }
        idx = by_site.get((contig, position))
        if idx is None:
            # A candidate with no cache row anywhere is a contract violation:
            # the candidate list was derived from this very cache.
            raise ValueError(
                f"{mag_id} {contig}:{position}: candidate rise but no cache row at the site"
            )
        site = rows.iloc[idx]
        site = pd.concat(
            [site, assess_rows(site, allele, model, min_freq, min_cov)], axis=1
        )
        later = site[
            (site["time"] == tps.later)
            & (site["group"] == group)
            & (site["replicate"] == replicate)
        ]
        carriers = sorted(set(later.loc[later["allele_present"], "subjectID"]))
        rec["carriers"] = ",".join(carriers)
        rec["n_carriers"] = len(carriers)
        if not carriers:
            rec["origin"] = ORIGIN_NOT_PRESENT_LATER
            records.append(rec)
            continue
        # own_group judges on the cage's group only; the other scopes see every
        # loaded sample and baseline_rows picks the earlier ones they define.
        pool = site[site["group"] == group] if scope == "own_group" else site
        t0 = baseline_rows(pool, scope, tps, set(carriers), baseline_cutoff(later))
        t0_cov = t0[t0["allele_status"] != EVIDENCE_NOT_COVERED]
        rec["n_baseline_samples_covered"] = int(len(t0_cov))
        rec["n_baseline_samples_present"] = int(t0_cov["allele_present"].sum())
        if t0_cov.empty:
            rec["origin"] = ORIGIN_BASELINE_NOT_COVERED
        elif t0_cov["allele_present"].any():
            rec["origin"] = ORIGIN_STANDING
        elif (t0_cov["allele_status"] == EVIDENCE_BELOW_DETECTION).any():
            rec["origin"] = ORIGIN_DE_NOVO_BELOW_DETECTION
        else:
            rec["origin"] = ORIGIN_DE_NOVO
        records.append(rec)

    out = pd.concat([rises.reset_index(drop=True), pd.DataFrame(records)], axis=1)

    # --- stage 4: blocks within a cage.  assign_blocks clusters by rose_in,
    # so the cage key keeps two cages' rises apart even at the same positions.
    alleles = out[["contig", "position", "allele", "origin"]].copy()
    alleles["mag_id"] = mag_id
    alleles["rose_in"] = [
        _cage_key(g, r)
        for g, r in zip(out["group"].astype(str), out["replicate"].astype(str))
    ]
    blocked = assign_blocks(alleles, window_bp, max_neighbours)
    out["origin"] = blocked["origin"].to_numpy()
    out["block_id"] = blocked["block_id"].to_numpy()
    out["n_de_novo_neighbours"] = blocked["n_de_novo_neighbours"].to_numpy()

    # --- stage 5: the carrier gate, de novo rows only, after the blocks.
    if strain_gate == "carrier":
        period = f"{tps.earlier}_{tps.later}"
        gated = out.index[out["origin"].isin(DE_NOVO_ORIGINS)]
        for i in gated:
            mice = [m for m in str(out.at[i, "carriers"]).split(",") if m]
            if turnover is None:
                status, n_det, n_chg = CARRIER_UNDETERMINED, 0, 0
            else:
                status, n_det, n_chg = carrier_strain_verdicts(
                    mice, turnover, period, min_determined_carriers
                )
            out.at[i, "carrier_strain_status"] = status
            out.at[i, "n_carriers_determined"] = n_det
            out.at[i, "n_carriers_strain_changed"] = n_chg
            if status == CARRIER_NEW_STRAIN:
                out.at[i, "origin"] = ORIGIN_NEW_STRAIN
            elif status == CARRIER_UNDETERMINED:
                out.at[i, "origin"] = CARRIER_UNDETERMINED
    return out[out_columns]


def is_event(origin: pd.Series, block_rule: str) -> pd.Series:
    """Which origins count as a de novo event for the gene tables.

    ``de_novo`` and ``de_novo_below_detection`` always; ``block_event`` only
    under ``block_rule == "once"`` (one event per block; ``exclude``, the
    default, drops blocks as Phase 2 does).  Nothing else ever counts.

    Example: origins [de_novo, block_event, standing_variation] -> exclude:
    [True, False, False]; once: [True, True, False].
    """
    if block_rule not in BLOCK_RULES:
        raise ValueError(f"block_rule must be one of {BLOCK_RULES}, got {block_rule!r}")
    counted = set(DE_NOVO_ORIGINS) | ({BLOCK_EVENT} if block_rule == "once" else set())
    return origin.isin(counted)


# ---------------------------------------------------------------------------
# Task 3: gene tables
# ---------------------------------------------------------------------------

SUBSTITUTION_CLASSES = ("synonymous", "nonsynonymous", "nonsense")
# The same complement table dN/dS uses for reverse-strand genes, so the two
# modules cannot drift apart on strand handling.
_COMPLEMENT = str.maketrans(
    "".join(COMPLEMENT_MAP), "".join(COMPLEMENT_MAP[b] for b in COMPLEMENT_MAP)
)
# Prodigal calls genes with the bacterial / archaeal table, as dN/dS does.
TRANSLATION_TABLE = 11

GENE_CAGE_COLUMNS = [
    "mag_id",
    "gene_id",
    "contig",
    "group",
    "replicate",
    "n_de_novo",
    "n_block_excluded",
    "n_synonymous",
    "n_nonsynonymous",
    "n_nonsense",
    "covered_bp",
    "gene_length",
    "covered_fraction",
    "entered",
    "rate",
]
CAGE_TOTAL_COLUMNS = [
    "mag_id",
    "group",
    "replicate",
    "n_de_novo_total",
    "covered_bp_total",
]


def classify_substitution(position: int, ref: str, alt: str, gene: pd.Series) -> str:
    """Does swapping ``ref`` for ``alt`` at this contig position change the protein?

    ``gene`` is a ``load_genes`` row: ``start, end`` (0-based inclusive),
    ``strand`` (1 / -1) and ``sequence`` (coding strand, 5' to 3').  The
    coding index comes from the shared ``coding_index`` (the same arithmetic
    dN/dS's ``get_codon_from_site`` routes through): ``position - start`` on
    the forward strand, ``end - position`` on the reverse strand, where the
    contig letters are complemented first (``COMPLEMENT_MAP``, shared with
    dN/dS) because the coding sequence is the reverse complement of the
    contig.  The codon is translated with the bacterial table before and
    after the swap.

    Raises when the position is outside the gene (named in the message), the
    contig's reference letter is not the coding sequence's letter there (the
    gene file and the reference disagree, or the strand is wrong), or the
    codon is incomplete (a sequence length that is not a multiple of three).

    Example: gene ATG GGA TGG at 100..108 forward; position 105 is the third
    base of GGA, A -> G gives GGG, still glycine: ``synonymous``; position 107
    G -> A turns TGG into TAG, a stop: ``nonsense``.
    """
    start, end, strand = int(gene["start"]), int(gene["end"]), int(gene["strand"])
    seq = str(gene["sequence"])
    # Re-raise the shared helper's out-of-gene error with the gene id in it:
    # callers see thousands of genes and need to know which one disagreed.
    if not start <= position <= end:
        raise ValueError(
            f"position {position} is outside gene {gene['gene_id']} ({start}..{end})"
        )
    index = coding_index(position, start, end, strand)
    if strand == 1:
        ref_c, alt_c = ref, alt
    else:
        ref_c, alt_c = ref.translate(_COMPLEMENT), alt.translate(_COMPLEMENT)
    if seq[index] != ref_c:
        raise ValueError(
            f"gene {gene['gene_id']} position {position}: the coding sequence has "
            f"{seq[index]!r} where the reference says {ref!r} (strand {strand})"
        )
    codon_start = (index // 3) * 3
    codon = seq[codon_start : codon_start + 3]
    if len(codon) != 3:
        raise ValueError(
            f"gene {gene['gene_id']} position {position}: partial codon {codon!r}"
        )
    offset = index - codon_start
    mutant = codon[:offset] + alt_c + codon[offset + 1 :]
    before = str(Seq(codon).translate(table=TRANSLATION_TABLE))
    after = str(Seq(mutant).translate(table=TRANSLATION_TABLE))
    if after == before:
        return "synonymous"
    if after == "*":
        return "nonsense"
    return "nonsynonymous"


def gene_tables(
    events: pd.DataFrame,
    readable: dict[tuple[str, str], dict[str, np.ndarray]],
    genes: pd.DataFrame,
    cages: pd.DataFrame,
    *,
    mag_id: str,
    block_rule: str,
    min_gene_covered_fraction: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The two tables the tests read: one row per gene x cage, one row per cage.

    ``events`` is Task 2's classified rise table plus ``ref_base`` (the
    contig's reference letter at the site, from the cache); ``is_event`` with
    ``block_rule`` decides which rows count.  ``readable`` is
    ``CoverageMask.readable``; ``genes`` is ``load_genes``; ``cages`` lists
    the comparison's (group, replicate) pairs so a cage with nothing readable
    still gets its zero rows.

    Gene x cage: ``n_de_novo`` (counted events in the gene), ``n_block_excluded``
    (block rows dropped by ``block_rule: exclude``), the three substitution
    classes of the counted events, ``covered_bp`` (readable bases inside the
    gene; the denominator, decision 7), ``gene_length``, ``covered_fraction =
    covered_bp / gene_length``, ``entered`` (``covered_bp > 0`` and the
    fraction at least ``min_gene_covered_fraction``, decision 5) and ``rate =
    n_de_novo / covered_bp`` (NaN when not entered).  A cage that does not
    enter for a gene reports ``covered_bp 0`` and ``n_de_novo 0`` so the tests
    never see it, while ``covered_fraction`` keeps the real number.

    Cage: ``n_de_novo_total`` = every counted event of the cage, inside a gene
    or not (an intergenic event is part of the genome's rate), and
    ``covered_bp_total``.

    An event whose ``gene_id`` is not in ``genes`` raises.

    Worked example (the unit test ``TestGeneTables``).  One 1,200 bp contig
    with three genes: g1 100..399 (+), g2 400..699 (-), g3 700..999 (+).
    Three cages: treatment A and B, control C.  Readable: A all 1,200 bases,
    B bases 0..599 only, C all 1,200.  Counted events: A at 150 and 160 (both
    in g1; A->G at a wobble position = synonymous, T->C at a first codon base
    = nonsynonymous) and at 1,100 (between genes); B at 450 (in g2, on the
    reverse strand; nonsynonymous).  C caught nothing.

    Gene x cage rows (9 in all; the zero rows are kept):

        gene  cage  n_de_novo  n_syn  n_nonsyn  covered_bp  fraction  entered  rate
        g1    A     2          1      1         300         1.00      True     2/300
        g2    A     0          0      0         300         1.00      True     0
        g2    B     1          0      1         200         0.67      True     1/200
        g3    B     0          0      0         0           0.00      False    NaN
        g1    C     0          0      0         300         1.00      True     0

    With ``min_gene_covered_fraction`` 0.8 the g2/B row no longer enters:
    ``entered`` False, ``covered_bp`` 0, ``n_de_novo`` 0, rate NaN, while
    ``covered_fraction`` still reads 0.67.

    Cage rows: A 3 events over 1,200 bases (the intergenic one counts), B 1
    over 600, C 0 over 1,200.  These are Test A's "fair share" numbers: cage
    A caught 3 in 1,200 bases, so a 300-base gene deserves 0.75 by chance.
    """
    if not 0.0 <= min_gene_covered_fraction <= 1.0:
        raise ValueError(
            f"min_gene_covered_fraction must be within [0, 1], got {min_gene_covered_fraction}"
        )
    gene_index = genes.set_index(genes["gene_id"].astype(str))
    cage_keys = [(str(g), str(r)) for g, r in zip(cages["group"], cages["replicate"])]

    # --- readable bases per gene and per cage; cages absent from the mask
    # (no paired sample at all) get zeros.
    known = {k: readable[k] for k in cage_keys if k in readable}
    by_gene = (
        covered_bp_by_gene(known, genes).set_index(["group", "replicate", "gene_id"])[
            "covered_bp"
        ]
        if known
        else pd.Series(dtype="int64")
    )
    totals_bp = (
        covered_bp_total(known).set_index(["group", "replicate"])["covered_bp_total"]
        if known
        else pd.Series(dtype="int64")
    )

    # --- the counted events, classified by substitution where they sit in a gene.
    ev = (
        events.copy()
        if len(events)
        else pd.DataFrame(
            columns=[
                "origin",
                "gene_id",
                "group",
                "replicate",
                "position",
                "ref_base",
                "allele",
            ]
        )
    )
    for col in ("gene_id", "group", "replicate", "ref_base", "allele"):
        ev[col] = ev[col].fillna("").astype(str) if col in ev.columns else ""
    counted = ev[is_event(ev["origin"], block_rule)].copy()
    excluded_blocks = ev[
        (ev["origin"] == BLOCK_EVENT) & ~is_event(ev["origin"], block_rule)
    ]
    classes = []
    for e in counted.itertuples(index=False):
        if e.gene_id == "":
            classes.append("")
            continue
        if e.gene_id not in gene_index.index:
            raise ValueError(
                f"event at {e.position} names gene {e.gene_id!r}, not in the gene file"
            )
        classes.append(
            classify_substitution(
                int(e.position), e.ref_base, e.allele, gene_index.loc[e.gene_id]
            )
        )
    counted["substitution"] = classes

    # --- gene x cage rows.
    records = []
    for group, replicate in cage_keys:
        in_cage = counted[
            (counted["group"] == group) & (counted["replicate"] == replicate)
        ]
        blocks_in_cage = excluded_blocks[
            (excluded_blocks["group"] == group)
            & (excluded_blocks["replicate"] == replicate)
        ]
        for gene in genes.itertuples(index=False):
            gid = str(gene.gene_id)
            length = int(gene.end) - int(gene.start) + 1
            covered = int(by_gene.get((group, replicate, gid), 0))
            fraction = covered / length if length else 0.0
            entered = covered > 0 and fraction >= min_gene_covered_fraction
            mine = in_cage[in_cage["gene_id"] == gid]
            n = int(len(mine)) if entered else 0
            sub = mine["substitution"] if entered else pd.Series(dtype=object)
            records.append(
                {
                    "mag_id": mag_id,
                    "gene_id": gid,
                    "contig": str(gene.contig),
                    "group": group,
                    "replicate": replicate,
                    "n_de_novo": n,
                    "n_block_excluded": int((blocks_in_cage["gene_id"] == gid).sum()),
                    "n_synonymous": int((sub == "synonymous").sum()),
                    "n_nonsynonymous": int((sub == "nonsynonymous").sum()),
                    "n_nonsense": int((sub == "nonsense").sum()),
                    "covered_bp": covered if entered else 0,
                    "gene_length": length,
                    "covered_fraction": fraction,
                    "entered": bool(entered),
                    "rate": (n / covered) if entered else np.nan,
                }
            )
    gene_cage = pd.DataFrame(records, columns=GENE_CAGE_COLUMNS)

    # --- cage totals: every counted event (intergenic included) over every readable base.
    totals = pd.DataFrame(
        [
            {
                "mag_id": mag_id,
                "group": group,
                "replicate": replicate,
                "n_de_novo_total": int(
                    (
                        (counted["group"] == group)
                        & (counted["replicate"] == replicate)
                    ).sum()
                ),
                "covered_bp_total": int(totals_bp.get((group, replicate), 0)),
            }
            for group, replicate in cage_keys
        ],
        columns=CAGE_TOTAL_COLUMNS,
    )
    return gene_cage, totals
