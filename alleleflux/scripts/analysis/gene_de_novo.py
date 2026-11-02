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

from alleleflux.scripts.analysis.allele_frequency._allele_freq_common import (
    NUCLEOTIDES,
)

logger = logging.getLogger(__name__)

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
