"""One test family -> its cells -> one MAG x cell grid (the pure half of the score figures).

A *family* is one statistic drawn consistently across the cells of a comparison:
parallelism in group A, parallelism in group B, divergence between them.  Six
families: ``ttest`` / ``wilcoxon`` (paired divergence), ``ttest_unpaired`` /
``mannwhitney`` (unpaired divergence, sharing the within-group cells with their
paired siblings), ``lmm`` and ``cmh``.  The
pipeline writes each cell to its own combined score table with its own column
name; ``FAMILY_SPECS`` is the single place that knows which file and which
column each cell lives in, and which ``(test_type, group_analyzed)`` row of the
significant-sites summary carries that cell's smallest q.

Nothing here draws or writes.  ``load_family`` returns two aligned frames,
``scores`` and ``min_q`` (rows = MAGs, columns = cell labels, NaN = untested),
plus the MAGs' taxonomy; the drawing module turns them into figures.

Example (family ``ttest``, comparison ``pre_end-treatment_control``):
``scores_single_sample-pre_end-treatment_control_treatment-MAGs.tsv`` column
``score_tTest_treatment (%)`` is the "treatment parallelism" cell, keyed
``("single_sample_tTest", "treatment")`` in the summary;
``scores_two_sample_paired-pre_end-treatment_control-MAGs.tsv`` column
``score_tTest (%)`` is the "divergence" cell, keyed ``("two_sample_paired_tTest", "")``.
"""

from __future__ import annotations

import glob
import logging
import os
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

TAXONOMY_COLUMNS = ["domain", "phylum", "class", "order", "family", "genus", "species"]

# family -> where each cell lives.  ``{group}`` is filled per parallelism group;
# the between-group (divergence) file has no group in its name.  CMH's
# divergence table carries the focus timepoint in its name (one file per focus),
# and its summary rows are keyed by that timepoint (``group_analyzed``).
FAMILY_SPECS = {
    "ttest": {
        "within_stem": "single_sample",
        "within_score": "score_tTest_{group} (%)",
        "within_test_type": "single_sample_tTest",
        "between_stem": "two_sample_paired",
        "between_score": "score_tTest (%)",
        "between_test_type": "two_sample_paired_tTest",
    },
    # the unpaired between-group tests share the within-group cells with their
    # paired siblings; they are the divergence test when samples cannot be paired
    # (subjects differ between timepoints, or a single-timepoint run)
    "ttest_unpaired": {
        "within_stem": "single_sample",
        "within_score": "score_tTest_{group} (%)",
        "within_test_type": "single_sample_tTest",
        "between_stem": "two_sample_unpaired",
        "between_score": "score_tTest (%)",
        "between_test_type": "two_sample_unpaired_tTest",
    },
    "wilcoxon": {
        "within_stem": "single_sample",
        "within_score": "score_Wilcoxon_{group} (%)",
        "within_test_type": "single_sample_Wilcoxon",
        "between_stem": "two_sample_paired",
        "between_score": "score_Wilcoxon (%)",
        "between_test_type": "two_sample_paired_Wilcoxon",
    },
    "mannwhitney": {
        "within_stem": "single_sample",
        "within_score": "score_Wilcoxon_{group} (%)",
        "within_test_type": "single_sample_Wilcoxon",
        "between_stem": "two_sample_unpaired",
        "between_score": "score_MannWhitney (%)",
        "between_test_type": "two_sample_unpaired_MannWhitney",
    },
    "lmm": {
        "within_stem": "lmm_across_time",
        "within_score": "score_LMM (%)",
        "within_test_type": "LMM_across_time",
        "between_stem": "lmm",
        "between_score": "score_LMM (%)",
        "between_test_type": "LMM",
    },
    "cmh": {
        "within_stem": "cmh_across_time",
        "within_score": "score_CMH (%)",
        "within_test_type": "CMH_across_time",
        "between_stem": "cmh",
        "between_score": "score_CMH (%)",
        "between_test_type": "CMH",
    },
}
FAMILIES = tuple(FAMILY_SPECS)


@dataclass(frozen=True)
class Cell:
    """One column of a family's grid.

    ``score_file`` is None when the pipeline wrote no table for the cell (no MAG
    was eligible): the cell is drawn as an all-grey column.  ``group_analyzed``
    is ``""`` for a between-group cell and the focus timepoint for CMH divergence,
    matching the summary table's key.
    """

    label: str  # column label, e.g. "treatment parallelism", "divergence (end)"
    kind: str  # "parallelism" | "divergence"
    group: str | None  # the parallelism group; None for divergence
    score_file: str | None
    score_column: str
    test_type: str
    group_analyzed: str


@dataclass(frozen=True)
class Family:
    name: str
    cells: tuple[Cell, ...]


@dataclass
class FamilyData:
    """``scores`` and ``min_q``: index = MAG id, columns = cell labels, NaN = untested.
    ``taxonomy``: index = MAG id, columns = TAXONOMY_COLUMNS."""

    family: Family
    scores: pd.DataFrame
    min_q: pd.DataFrame
    taxonomy: pd.DataFrame


def total_column(score_column: str) -> str:
    """``"score_X (%)"`` -> ``"total_sites_per_group_X"`` (the sites behind the score).

    Example: ``score_tTest_control (%)`` -> ``total_sites_per_group_tTest_control``.
    """
    core = score_column[len("score_") : -len(" (%)")]
    return f"total_sites_per_group_{core}"


def resolve_family(
    scores_dir: str,
    comparison: str,
    groups: Sequence[str],
    family: str,
    data_type: str = "longitudinal",
) -> Family | None:
    """The cells of ``family`` for one comparison, with the files that exist.

    Parallelism cells, one per group in ``groups`` (longitudinal only: a single
    timepoint has no within-group test), then the divergence cell(s).  A cell
    whose file is missing keeps ``score_file=None``.  Returns None when no cell
    has a file, so the caller can skip the family with one log line.

    Example: groups ["treatment", "control"], family "ttest" -> three cells
    labelled "treatment parallelism", "control parallelism", "divergence".
    """
    spec = FAMILY_SPECS[family]
    cells: list[Cell] = []

    def _existing(path: str) -> str | None:
        return path if os.path.exists(path) else None

    if data_type == "longitudinal":
        for group in groups:
            path = os.path.join(
                scores_dir,
                f"scores_{spec['within_stem']}-{comparison}_{group}-MAGs.tsv",
            )
            cells.append(
                Cell(
                    label=f"{group} parallelism",
                    kind="parallelism",
                    group=group,
                    score_file=_existing(path),
                    score_column=spec["within_score"].format(group=group),
                    test_type=spec["within_test_type"],
                    group_analyzed=group,
                )
            )
    if family == "cmh":
        # one divergence table per focus timepoint, keyed by that timepoint
        pattern = os.path.join(scores_dir, f"scores_cmh-{comparison}-MAGs-*.tsv")
        focus_files = sorted(glob.glob(pattern))
        if not focus_files:
            cells.append(
                Cell(
                    "divergence",
                    "divergence",
                    None,
                    None,
                    spec["between_score"],
                    "CMH",
                    "",
                )
            )
        for path in focus_files:
            focus = path[: -len(".tsv")].rsplit("-MAGs-", 1)[1]
            cells.append(
                Cell(
                    label=f"divergence ({focus})",
                    kind="divergence",
                    group=None,
                    score_file=path,
                    score_column=spec["between_score"],
                    test_type="CMH",
                    group_analyzed=focus,
                )
            )
    else:
        path = os.path.join(
            scores_dir, f"scores_{spec['between_stem']}-{comparison}-MAGs.tsv"
        )
        cells.append(
            Cell(
                label="divergence",
                kind="divergence",
                group=None,
                score_file=_existing(path),
                score_column=spec["between_score"],
                test_type=spec["between_test_type"],
                group_analyzed="",
            )
        )
    if not any(c.score_file for c in cells):
        return None
    return Family(family, tuple(cells))


def load_family(
    family: Family, cell_stats: pd.DataFrame, comparison: str
) -> FamilyData:
    """Read every cell's table into one grid; look up each cell's smallest q.

    A MAG is *tested* in a cell when its row exists and the sites behind the
    score are not zero (the CMH combiner keeps zero-site MAGs as score 0,
    "insufficient data"; those are untested here).  Rows = the union of tested
    MAGs over the cells; an untested cell is NaN in both frames.  ``min_q`` comes
    from ``cell_stats`` rows matching ``(comparison, test_type, group_analyzed)``.
    Raises ValueError when a score table lacks its total-sites column, when a
    cell has duplicate summary rows, or when the tested MAGs of a cell differ
    between the score table and the summary (both derive from the same p-value
    file, so any difference means a broken run).

    Example: M2 tested only in "treatment parallelism" -> row M2 with one score
    and NaN in the other columns, in both frames.
    """
    scores, tax_frames = {}, []
    for cell in family.cells:
        if cell.score_file is None:
            scores[cell.label] = pd.Series(dtype=float)
            continue
        df = pd.read_csv(cell.score_file, sep="\t", dtype={"MAG_ID": str})
        total = total_column(cell.score_column)
        if total not in df.columns:
            # without the site count we cannot tell a real 0 from "not tested"
            raise ValueError(f"{cell.score_file}: missing column {total}")
        df = df[df[total].astype(float) > 0]
        scores[cell.label] = df.set_index("MAG_ID")[cell.score_column].astype(float)
        present = [c for c in TAXONOMY_COLUMNS if c in df.columns]
        tax_frames.append(df.set_index("MAG_ID")[present])
    score_df = pd.DataFrame(scores)  # outer join on MAG id; NaN where untested
    score_df.index.name = "MAG_ID"
    score_df = score_df[[c.label for c in family.cells]]

    stats = cell_stats[cell_stats["comparison"] == comparison]
    key_group = stats["group_analyzed"].fillna("").astype(str)
    min_q = pd.DataFrame(index=score_df.index, columns=score_df.columns, dtype=float)
    for cell in family.cells:
        rows = stats[
            (stats["test_type"] == cell.test_type) & (key_group == cell.group_analyzed)
        ]
        q = rows.set_index("mag_id")["min_q_value"].astype(float)
        dup = q.index[q.index.duplicated()]
        if len(dup):
            raise ValueError(
                f"{family.name} / {cell.label}: duplicate summary rows for {list(dup)[:5]}"
            )
        # the score table and the summary are built from the same per-MAG p-value
        # file, so the tested MAG set must be identical on both sides
        tested = set(score_df.index[score_df[cell.label].notna()])
        no_q = sorted(tested - set(q.index))
        no_score = sorted(set(q.index) - tested)
        if no_q or no_score:
            raise ValueError(
                f"{family.name} / {cell.label}: score/summary mismatch; "
                f"scored but no summary row: {no_q[:5]}; "
                f"summary row but no score: {no_score[:5]}"
            )
        min_q[cell.label] = q.reindex(score_df.index)

    taxonomy = (
        pd.concat(tax_frames).groupby(level=0).first().reindex(score_df.index)
        if tax_frames
        else pd.DataFrame(index=score_df.index)
    )
    taxonomy = taxonomy.reindex(columns=TAXONOMY_COLUMNS)
    return FamilyData(family, score_df, min_q, taxonomy)


def taxonomy_order(taxonomy: pd.DataFrame) -> list[str]:
    """MAG ids sorted by phylum, class, order, family, genus, species, then id; blanks last."""
    keys = [c for c in TAXONOMY_COLUMNS[1:] if c in taxonomy.columns]
    ordered = taxonomy.assign(_id=taxonomy.index.astype(str)).sort_values(
        keys + ["_id"], na_position="last"
    )
    return list(ordered.index)


def tree_order(mags: Iterable[str], newick_path: str) -> tuple[list[str], list[str]]:
    """MAG ids in the tree's (ladderized) tip order, and those absent from the tree.

    Example: tree ``((M2,M1),M3);`` and mags [M1, M2, M3, M4] -> ([M2, M1, M3]
    up to ladderizing, [M4]).
    """
    from Bio import Phylo

    tree = Phylo.read(newick_path, "newick")
    tree.ladderize()
    tips = [t.name for t in tree.get_terminals()]
    wanted = list(mags)
    in_tree = [t for t in tips if t in set(wanted)]
    absent = [m for m in wanted if m not in set(tips)]
    return in_tree, absent


def q_stars(q: float, threshold: float) -> str:
    """Significance glyph for a q-value: *** < 0.001, ** < 0.01, * < threshold, else ''."""
    if q is None or np.isnan(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < threshold:
        return "*"
    return ""
