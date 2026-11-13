"""Gene-level de novo accumulation: which genes collect new alleles, per cage.

The question (spec section 3).  Every cage keeps a diary of the new letters
that appeared in a genome between the two timepoints of a comparison.  Each
entry is filed under the gene it landed in.  Test A asks whether a gene got
more entries than its share of the readable genome inside one group; Test B
asks whether a gene collected entries per readable base faster in one group
than in the other.  The cage is the replicate throughout.

This module holds the COUNTING side as pure functions over frames and the
command line that streams one comparison's files through them
(``alleleflux-gene-de-novo``: ``hunt_the_hot_genes`` orchestrates, one
``_count_one_mag`` worker per MAG); the two Poisson tests live in
``statistics/poisson_rate_tests.py``.  Four files come out per comparison:
``{mag}_gene_de_novo_events.tsv`` (every candidate rise with its verdict),
``{mag}_gene_de_novo_genes.tsv`` (gene x cage), ``{mag}_gene_de_novo_cages.tsv``
(cage totals) and ``gene_de_novo_tests_{comparison}.tsv`` (one row per MAG x
gene, both tests, BH q-values).  The counting functions:

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

import argparse
import glob
import logging
import multiprocessing
import os

import numpy as np
import pandas as pd

from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
    NUCLEOTIDES,
)
from alleleflux.scripts.analysis.ani.null_model import build_error_model
from alleleflux.scripts.analysis.baseline_presence import (
    EVIDENCE_BELOW_DETECTION,
    EVIDENCE_NOT_COVERED,
    Timepoints,
    assess_rows,
    baseline_cutoff,
    baseline_rows,
    parse_comparison,
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
from alleleflux.scripts.statistics.poisson_rate_tests import (
    result_columns,
    test_every_gene,
)
from alleleflux.scripts.utilities.cache_io import iter_cache_batches, load_cache_sites
from alleleflux.scripts.utilities.genes import load_genes  # noqa: F401 (re-export)
from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.utilities.sample_metadata import (
    META_OPTIONAL,
    validate_collection_dates,
)
from alleleflux.scripts.utilities.utilities import build_contig_length_index

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
# The events file: every candidate rise with its verdict (``classify_cage_rises``).
EVENT_COLUMNS = RISE_COLUMNS + VERDICT_COLUMNS


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


def _cohort_cutoff(
    dates: np.ndarray, cage_later: np.ndarray, status: np.ndarray, present: np.ndarray
):
    """``baseline_presence.baseline_cutoff`` on arrays: the date one cage's baseline is cut at.

    ``dates`` are the site's rows' collection dates (NaT when blank),
    ``cage_later`` flags the cage's later-timepoint rows, ``status`` /
    ``present`` are the letter's verdicts.  The first carrier's date when the
    cage has a covered carrier (None if any carrier is undated); otherwise the
    cage's earliest later sampling day (None if any of those rows is undated
    or there is none).  Returned as ``datetime64`` for a direct ``<``.

    Example: carriers dated 2024-04-15 and 2024-05-01 -> 2024-04-15; the
    second one blank -> None.
    """
    t1_cov = cage_later & (status != EVIDENCE_NOT_COVERED)
    carriers = dates[t1_cov & present]
    if len(carriers):
        return None if np.isnat(carriers).any() else carriers.min()
    later = dates[cage_later]
    if not len(later) or np.isnat(later).any():
        return None
    return later.min()


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

    # --- the rows every verdict is read from, as flat arrays.  The presence
    # rule (assess_rows, already vectorised over rows) runs ONCE per letter over
    # every loaded row instead of once per candidate: 70 k candidates through
    # per-site frames took 167 s, the same verdicts in bulk take about a second.
    # Labels become text (the cache stores categoricals), dates come from the
    # sheet by sample id (the cohort baseline is cut by date), and each site's
    # row positions are indexed once so a candidate is a few array lookups.
    rows = cache_rows.reset_index(drop=True)
    contig = rows["contig"].astype(str).to_numpy()
    position = rows["position"].astype("int64").to_numpy()
    subject = rows["subjectID"].astype(str).to_numpy()
    replicate = rows["replicate"].astype(str).to_numpy()
    group = rows["group"].astype(str).to_numpy()
    time = rows["timepoint"].astype(str).to_numpy()
    meta = sheet.set_index(sheet["sample_id"].astype(str))
    dated = (
        pd.to_datetime(meta["collection_date"], format="%Y-%m-%d", errors="coerce")
        if "collection_date" in meta.columns
        else pd.Series(pd.NaT, index=meta.index)
    )
    dates = dated.reindex(rows["sample_id"].astype(str)).to_numpy(
        dtype="datetime64[ns]"
    )
    is_later = time == str(tps.later)
    is_earlier = time == str(tps.earlier)
    # One verdict column per letter that rose anywhere: (status words, present flags).
    verdict = {}
    for letter in sorted(set(rises["allele"].astype(str))):
        v = assess_rows(rows, letter, model, min_freq, min_cov)
        verdict[letter] = (
            v["allele_status"].to_numpy(),
            v["allele_present"].to_numpy(dtype=bool),
        )
    by_site = {
        key: np.asarray(idx)
        for key, idx in pd.DataFrame({"contig": contig, "position": position})
        .groupby(["contig", "position"])
        .indices.items()
    }

    # --- stages 1-3, one candidate at a time over the arrays.  The baseline
    # selection mirrors baseline_presence.baseline_rows / baseline_cutoff word
    # for word (the unit test TestClassifyInBulk pins the two against each
    # other on random sites), so the cage-level and site-level analyses can
    # never disagree on what "absent at baseline" means.
    records = []
    for rise in rises.itertuples(index=False):
        key = (str(rise.contig), int(rise.position))
        g, r, letter = str(rise.group), str(rise.replicate), str(rise.allele)
        rec = {
            "carriers": "",
            "n_carriers": 0,
            "carrier_strain_status": "",
            "n_carriers_determined": 0,
            "n_carriers_strain_changed": 0,
            "n_baseline_samples_covered": 0,
            "n_baseline_samples_present": 0,
        }
        idx = by_site.get(key)
        if idx is None:
            # A candidate with no cache row anywhere is a contract violation:
            # the candidate list was derived from this very cache.
            raise ValueError(
                f"{mag_id} {key[0]}:{key[1]}: candidate rise but no cache row at the site"
            )
        status, present = verdict[letter]
        s_idx, p_idx = status[idx], present[idx]
        # stage 2: the cage's later-timepoint mice that carry the letter.
        cage_later = is_later[idx] & (group[idx] == g) & (replicate[idx] == r)
        carriers = np.unique(subject[idx][cage_later & p_idx])
        rec["carriers"] = ",".join(carriers.tolist())
        rec["n_carriers"] = int(len(carriers))
        if not len(carriers):
            rec["origin"] = ORIGIN_NOT_PRESENT_LATER
            records.append(rec)
            continue
        # stage 3: which of the site's rows are the baseline (baseline_rows):
        if scope == "individual":
            base = is_earlier[idx] & np.isin(subject[idx], carriers)
        elif scope == "own_group":
            base = is_earlier[idx] & (group[idx] == g)
        elif scope == "comparison":
            base = is_earlier[idx]
        else:  # cohort: every row dated before the cutoff (baseline_cutoff);
            # an undated sheet narrows to the earlier timepoint, as baseline_rows does.
            d_idx = dates[idx]
            cutoff = _cohort_cutoff(d_idx, cage_later, s_idx, p_idx)
            if cutoff is not None and not np.isnat(d_idx).any():
                base = d_idx < cutoff
            else:
                base = is_earlier[idx]
        covered = base & (s_idx != EVIDENCE_NOT_COVERED)
        rec["n_baseline_samples_covered"] = int(covered.sum())
        rec["n_baseline_samples_present"] = int((covered & p_idx).sum())
        if not covered.any():
            rec["origin"] = ORIGIN_BASELINE_NOT_COVERED
        elif (covered & p_idx).any():
            rec["origin"] = ORIGIN_STANDING
        elif (covered & (s_idx == EVIDENCE_BELOW_DETECTION)).any():
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

GENE_CAGE_COLUMNS = [
    "mag_id",
    "gene_id",
    "contig",
    "group",
    "replicate",
    "n_de_novo",
    "n_block_excluded",
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

    ``events`` is Task 2's classified rise table; ``is_event`` with
    ``block_rule`` decides which rows count.  ``readable`` is
    ``CoverageMask.readable``; ``genes`` is ``load_genes``; ``cages`` lists
    the comparison's (group, replicate) pairs so a cage with nothing readable
    still gets its zero rows.

    Gene x cage: ``n_de_novo`` (counted events in the gene), ``n_block_excluded``
    (block rows dropped by ``block_rule: exclude``), ``covered_bp`` (readable bases inside the
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

    No synonymous / nonsynonymous split (decision 4, deferred 2026-11-03): a
    class computed against the assembly's codon is wrong wherever the
    population differs from the reference, which with an unrelated catalog
    genome is at every fixed difference.  When the split returns it must be
    computed against the baseline codon, the way dN/dS reconstructs its
    ancestral sequence, as extra columns that leave ``n_de_novo`` untouched.

    Worked example (the unit test ``TestGeneTables``).  One 1,200 bp contig
    with three genes: g1 100..399 (+), g2 400..699 (-), g3 700..999 (+).
    Three cages: treatment A and B, control C.  Readable: A all 1,200 bases,
    B bases 0..599 only, C all 1,200.  Counted events: A at 150 and 160 (both
    in g1) and at 1,100 (between genes); B at 450 (in g2).  C caught nothing.

    Gene x cage rows (9 in all; the zero rows are kept):

        gene  cage  n_de_novo  covered_bp  fraction  entered  rate
        g1    A     2          300         1.00      True     2/300
        g2    A     0          300         1.00      True     0
        g2    B     1          200         0.67      True     1/200
        g3    B     0          0           0.00      False    NaN
        g1    C     0          300         1.00      True     0

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

    # --- the counted events, each placed in its gene (or none).
    ev = (
        events.copy()
        if len(events)
        else pd.DataFrame(
            columns=["origin", "gene_id", "group", "replicate", "position"]
        )
    )
    for col in ("gene_id", "group", "replicate"):
        ev[col] = ev[col].fillna("").astype(str) if col in ev.columns else ""
    counted = ev[is_event(ev["origin"], block_rule)].copy()
    excluded_blocks = ev[
        (ev["origin"] == BLOCK_EVENT) & ~is_event(ev["origin"], block_rule)
    ]
    # An event naming a gene the gene file does not know is a contract
    # violation (the gene ids come from the same Prodigal file), never a
    # silent zero.
    named = counted[counted["gene_id"] != ""]
    unknown = named[~named["gene_id"].isin(gene_index.index)]
    if len(unknown):
        e = unknown.iloc[0]
        raise ValueError(
            f"event at {e.position} names gene {e.gene_id!r}, not in the gene file"
        )

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
            records.append(
                {
                    "mag_id": mag_id,
                    "gene_id": gid,
                    "contig": str(gene.contig),
                    "group": group,
                    "replicate": replicate,
                    "n_de_novo": n,
                    "n_block_excluded": int((blocks_in_cage["gene_id"] == gid).sum()),
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


# ---------------------------------------------------------------------------
# Task 6: the command line (one comparison, every MAG of it)
# ---------------------------------------------------------------------------

MEAN_CHANGE_SUFFIX = "_allele_frequency_changes_mean.parquet"


def _mags_in(allele_analysis_dir: str) -> list[str]:
    """Every MAG with a mean-change table in the comparison's allele-analysis folder."""
    paths = glob.glob(os.path.join(allele_analysis_dir, f"*{MEAN_CHANGE_SUFFIX}"))
    mags = sorted(os.path.basename(p)[: -len(MEAN_CHANGE_SUFFIX)] for p in paths)
    if not mags:
        raise ValueError(
            f"no *{MEAN_CHANGE_SUFFIX} in {allele_analysis_dir}: nothing to count"
        )
    return mags


def _per_mag_paths(output_dir: str, mag_id: str) -> dict[str, str]:
    return {
        kind: os.path.join(output_dir, f"{mag_id}_gene_de_novo_{kind}.tsv")
        for kind in ("events", "genes", "cages")
    }


def _read_turnover(turnover_dir: str, mag_id: str) -> pd.DataFrame | None:
    """This MAG's per-mouse strain verdicts, or None (then every event is undetermined)."""
    path = os.path.join(turnover_dir, f"{mag_id}_strain_turnover.tsv")
    if not os.path.exists(path):
        logger.warning(
            f"{mag_id}: no turnover table at {path}; every de novo event of this "
            f"MAG is strain_undetermined and not counted"
        )
        return None
    return pd.read_csv(path, sep="\t", dtype=str)


def _count_one_mag(job: dict) -> tuple[str, pd.DataFrame, pd.DataFrame, dict]:
    """Worker: one MAG from its mean-change table and cache to its three files.

    ``job`` carries the MAG id, its genes and contig lengths, the paths, the
    comparison and every knob (see ``hunt_the_hot_genes``).  Steps: candidate
    rises -> readable bases per cage (a streamed pass over the paired samples'
    cache rows, every position) -> the candidates' cache rows for every loaded
    sample -> verdicts -> gene x cage and cage tables -> written.  Returns the
    two tables for the tests plus the origin tally for the log.
    """
    mag_id = job["mag_id"]
    tps: Timepoints = job["tps"]
    mean_changes = pd.read_parquet(job["mean_change_path"])
    rises = cage_rises(mean_changes, job["groups"])

    # Readable bases: one bit per position per cage and timepoint, from the
    # paired mice only (decision 7), streamed so a large cache never sits in
    # memory at once.
    mask = CoverageMask(job["contig_lengths"], job["min_cov"])
    for tp in (tps.earlier, tps.later):
        for batch in iter_cache_batches(
            job["cache_dir"],
            mag_id,
            tp,
            sample_ids=job["paired_sample_ids"],
            columns=["total_coverage", "group", "replicate"],
        ):
            mask.add(batch)
    readable = mask.readable(tps.earlier, tps.later)

    if rises.empty:
        events = pd.DataFrame(columns=EVENT_COLUMNS)
    else:
        sites = rises[["contig", "position"]].drop_duplicates()
        rows = load_cache_sites(
            job["cache_dir"],
            mag_id,
            job["timepoints"],
            sites=sites,
            sample_ids=job["loaded_sample_ids"],
            columns=[
                "A",
                "C",
                "G",
                "T",
                "total_coverage",
                "subjectID",
                "replicate",
                "group",
            ],
        )
        turnover = (
            _read_turnover(job["turnover_dir"], mag_id)
            if job["strain_gate"] == "carrier"
            else None
        )
        events = classify_cage_rises(
            rises,
            rows,
            job["sheet"],
            mag_id=mag_id,
            tps=tps,
            scope=job["baseline_scope"],
            model=job["model"],
            min_freq=job["min_freq"],
            min_cov=job["min_cov"],
            window_bp=job["cluster_window_bp"],
            max_neighbours=job["max_de_novo_neighbours"],
            strain_gate=job["strain_gate"],
            min_determined_carriers=job["min_determined_carriers"],
            turnover=turnover,
        )
    gene_cage, totals = gene_tables(
        events,
        readable,
        job["genes"],
        job["cages"],
        mag_id=mag_id,
        block_rule=job["block_rule"],
        min_gene_covered_fraction=job["min_gene_covered_fraction"],
    )
    paths = _per_mag_paths(job["output_dir"], mag_id)
    events[EVENT_COLUMNS].to_csv(paths["events"], sep="\t", index=False)
    gene_cage.to_csv(paths["genes"], sep="\t", index=False)
    totals.to_csv(paths["cages"], sep="\t", index=False)
    tally = events["origin"].value_counts().to_dict() if len(events) else {}
    return mag_id, gene_cage, totals, tally


def hunt_the_hot_genes(args: argparse.Namespace) -> int:
    """Orchestrator: one comparison -> per-MAG tables -> one tests table.

    1. The sample sheet: the comparison's samples (both groups, both
       timepoints) decide the cages and, through ``paired_samples``, the
       denominators; under the cohort scope every sample of the sheet is
       loaded as potential baseline, and its dates must be all-or-nothing
       (a partly dated sheet is an error, an undated one narrows the baseline
       to the earlier timepoint with a warning), as in baseline presence.
    2. Genes and contig lengths once (``load_genes``,
       ``build_contig_length_index``), split per MAG through the MAG mapping.
    3. One worker per MAG (``_count_one_mag``), the MAGs being every
       mean-change table in ``--allele_analysis_dir``.
    4. ``poisson_rate_tests.test_every_gene`` over every MAG's tables ->
       ``gene_de_novo_tests_{comparison}.tsv`` (one row per MAG x gene).
    """
    earlier, later, group_a, group_b = parse_comparison(args.comparison)
    tps = Timepoints(earlier, later)
    groups = (group_a, group_b)
    os.makedirs(args.output_dir, exist_ok=True)

    # ---- 1. the sheet
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    if "replicate" not in meta.columns:
        meta["replicate"] = meta["subjectID"]  # same default as mag_metadata.py
    for col in META_OPTIONAL:
        if col not in meta.columns:
            meta[col] = np.nan
    validate_collection_dates(meta)
    in_comparison = meta["group"].isin(groups) & meta["time"].isin((earlier, later))
    comparison_meta = meta[in_comparison]
    if comparison_meta.empty:
        raise ValueError(
            f"no samples in {args.metadata} for groups {groups} at {earlier}/{later}"
        )
    for g in groups:
        if not (comparison_meta["group"] == g).any():
            raise ValueError(
                f"group {g!r} has no sample at {earlier}/{later} in {args.metadata}"
            )
    paired = paired_samples(comparison_meta, earlier, later)
    if paired.empty:
        raise ValueError(
            f"no mouse of groups {groups} is sampled at both {earlier} and {later}: "
            f"nothing was searched for a rise (check subjectID spelling)"
        )
    cages = (
        comparison_meta[["group", "replicate"]]
        .drop_duplicates()
        .sort_values(["group", "replicate"])
        .reset_index(drop=True)
    )
    logger.info(
        f"{len(paired)} paired samples ({paired['subjectID'].nunique()} mice) over "
        f"{len(cages)} cages: "
        + ", ".join(f"{g}={int((cages['group'] == g).sum())}" for g in groups)
    )
    # Which samples and timepoints the verdicts read: the comparison's own,
    # plus every other sample of the sheet under a DATED cohort scope.
    loaded = comparison_meta
    timepoints = [earlier, later]
    if args.baseline_scope == "cohort":
        dated = meta["collection_date"].notna() & (
            meta["collection_date"].astype(str).str.strip() != ""
        )
        if dated.any() and not dated.all():
            missing = sorted(meta.loc[~dated, "sample_id"])
            raise ValueError(
                f"--baseline_scope cohort needs a collection_date (YYYY-MM-DD) for "
                f"every sample or for none; {len(missing)} of {len(meta)} are blank: "
                f"{missing}"
            )
        if dated.all():
            loaded = meta
            timepoints = list(dict.fromkeys(meta["time"]))
            logger.info(
                f"cohort scope: {len(meta) - len(comparison_meta)} samples outside "
                f"the comparison loaded as potential baseline (timepoints {timepoints})"
            )
        else:
            logger.warning(
                f"cohort scope: no collection_date in the metadata, so the baseline "
                f"is the comparison's earlier timepoint {earlier!r} only (as under "
                f"--baseline_scope comparison).  Add a collection_date column to "
                f"use the full cohort."
            )

    # ---- 2. genes and contig lengths, once, then per MAG
    genes = load_genes(args.prodigal_fasta)
    mapping = pd.read_csv(args.mag_mapping, sep="\t", dtype=str)
    contig_lengths = build_contig_length_index(args.fasta, args.mag_mapping)
    contigs_of = mapping.groupby("mag_id")["contig_id"].apply(list).to_dict()
    mags = _mags_in(args.allele_analysis_dir)
    logger.info(
        f"{len(mags)} MAGs with a mean-change table in {args.allele_analysis_dir}; "
        f"{len(genes)} genes in {args.prodigal_fasta}"
    )
    model = build_error_model(min_base_quality=args.min_base_quality, fdr=args.fdr)
    jobs = []
    for mag_id in mags:
        contigs = contigs_of.get(mag_id)
        if not contigs:
            raise ValueError(f"MAG {mag_id} has no contig in {args.mag_mapping}")
        missing = [c for c in contigs if c not in contig_lengths]
        if missing:
            raise ValueError(
                f"MAG {mag_id}: contig(s) {missing[:5]} not in {args.fasta}"
            )
        jobs.append(
            {
                "mag_id": mag_id,
                "mean_change_path": os.path.join(
                    args.allele_analysis_dir, f"{mag_id}{MEAN_CHANGE_SUFFIX}"
                ),
                "genes": genes[genes["contig"].isin(contigs)].reset_index(drop=True),
                "contig_lengths": {c: contig_lengths[c] for c in contigs},
                "cages": cages,
                "sheet": meta,
                "groups": groups,
                "tps": tps,
                "timepoints": timepoints,
                "paired_sample_ids": list(paired["sample_id"]),
                "loaded_sample_ids": list(loaded["sample_id"]),
                "cache_dir": args.cache_dir,
                "turnover_dir": args.turnover_dir,
                "output_dir": args.output_dir,
                "model": model,
                "min_cov": args.min_cov,
                "min_freq": args.min_freq,
                "baseline_scope": args.baseline_scope,
                "cluster_window_bp": args.cluster_window_bp,
                "max_de_novo_neighbours": args.max_de_novo_neighbours,
                "strain_gate": args.strain_gate,
                "min_determined_carriers": args.min_determined_carriers,
                "block_rule": args.block_rule,
                "min_gene_covered_fraction": args.min_gene_covered_fraction,
            }
        )

    # ---- 3. one worker per MAG
    gene_frames, cage_frames = [], []
    with multiprocessing.Pool(processes=args.cpus) as pool:
        for mag_id, gene_cage, totals, tally in pool.imap_unordered(
            _count_one_mag, jobs
        ):
            gene_frames.append(gene_cage)
            cage_frames.append(totals)
            counted = int(totals["n_de_novo_total"].sum())
            logger.info(
                f"{mag_id}: {sum(tally.values()):,} candidate rises -> {tally}; "
                f"{counted} events counted over {len(totals)} cages"
            )

    # ---- 4. the tests, over every MAG of the comparison
    gene_cage = pd.concat(gene_frames, ignore_index=True)
    totals = pd.concat(cage_frames, ignore_index=True)
    tests = test_every_gene(
        gene_cage,
        totals,
        groups,
        min_cages=args.min_cages,
        fdr_by_mag=args.fdr_group_by_mag_id,
    )
    path = os.path.join(args.output_dir, f"gene_de_novo_tests_{args.comparison}.tsv")
    tests.to_csv(path, sep="\t", index=False)
    n_between = int(tests["p_value_between"].notna().sum()) if len(tests) else 0
    logger.info(
        f"wrote {len(tests):,} MAG x gene rows to {path}; {n_between:,} genes tested "
        f"between groups at the floor of {args.min_cages} cages per group"
        + (
            f"; smallest p_value_between {tests['p_value_between'].min():.3g} "
            f"({tests.loc[tests['p_value_between'].idxmin(), 'gene_id']})"
            if n_between
            else ""
        )
    )
    return 0


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Gene-level de novo accumulation for one comparison: which genes "
        "collect new alleles faster than their genome (per group) or than the other "
        "group, with the cage as the replicate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--comparison",
        required=True,
        help="Comparison label, e.g. pre_post-fat_control",
    )
    parser.add_argument(
        "--cache_dir",
        required=True,
        help="The run's allele_freq_cache/ directory (QC-passing samples, site-filtered)",
    )
    parser.add_argument(
        "--allele_analysis_dir",
        required=True,
        help="The comparison's allele_analysis folder: one "
        "{mag}_allele_frequency_changes_mean.parquet per MAG = the MAGs analysed",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Sample metadata TSV (sample_id, subjectID, group, time[, replicate, collection_date])",
    )
    parser.add_argument(
        "--prodigal_fasta", required=True, help="Prodigal gene file (.fna)"
    )
    parser.add_argument(
        "--fasta", required=True, help="Reference FASTA (contig lengths)"
    )
    parser.add_argument(
        "--mag_mapping",
        required=True,
        help="MAG-to-contig mapping TSV (mag_id, contig_id)",
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--turnover_dir",
        default=None,
        help="alleleflux-strain-turnover outputs ({mag}_strain_turnover.tsv); the carrier gate reads them",
    )
    parser.add_argument(
        "--baseline_scope",
        choices=BASELINE_SCOPES,
        default="cohort",
        help="What 'absent at baseline' is judged against (as alleleflux-baseline-presence)",
    )
    parser.add_argument(
        "--cluster_window_bp",
        type=int,
        default=1000,
        help="cluster rule: de novo alleles within this many bp on one contig IN ONE CAGE are neighbours",
    )
    parser.add_argument(
        "--max_de_novo_neighbours",
        type=int,
        default=0,
        help="cluster rule: more de novo neighbours than this -> block_event",
    )
    parser.add_argument(
        "--strain_gate",
        choices=STRAIN_GATES,
        default="carrier",
        help="carrier = the event's own carriers' turnover verdicts (needs --turnover_dir); none = no gate",
    )
    parser.add_argument(
        "--min_determined_carriers",
        type=int,
        default=1,
        help="carrier gate: carriers with a verdict needed; fewer -> strain_undetermined",
    )
    parser.add_argument(
        "--block_rule",
        choices=BLOCK_RULES,
        default="exclude",
        help="block events: exclude (not counted) or once (one event per block row)",
    )
    parser.add_argument(
        "--min_gene_covered_fraction",
        type=float,
        default=0.0,
        help="a cage enters a gene only when at least this fraction of it is readable",
    )
    parser.add_argument(
        "--min_cages",
        type=int,
        required=True,
        help="cages per group a gene needs for a test (quality_control.min_sample_num)",
    )
    parser.add_argument(
        "--fdr_group_by_mag_id",
        action="store_true",
        help="Benjamini-Hochberg per MAG instead of over every gene of the comparison",
    )
    parser.add_argument(
        "--min_cov", type=int, default=5, help="Reads before any verdict"
    )
    parser.add_argument(
        "--min_freq", type=float, default=0.05, help="Presence frequency floor"
    )
    parser.add_argument(
        "--fdr", type=float, default=1e-6, help="Null-model false-discovery rate"
    )
    parser.add_argument(
        "--min_base_quality",
        type=int,
        default=30,
        help="Base quality the null model assumes",
    )
    parser.add_argument("--cpus", type=int, default=multiprocessing.cpu_count())
    args = parser.parse_args()
    if args.strain_gate == "carrier" and not args.turnover_dir:
        parser.error("--strain_gate carrier needs --turnover_dir")
    if args.min_cages < 1:
        parser.error("--min_cages must be >= 1")
    return hunt_the_hot_genes(args)


if __name__ == "__main__":
    raise SystemExit(main())
