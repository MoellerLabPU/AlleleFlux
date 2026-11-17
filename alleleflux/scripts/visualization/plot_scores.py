#!/usr/bin/env python3
"""Score heatmaps and boxplots for one comparison: ``alleleflux-plot-scores``.

Reads a finished run's combined MAG score tables and its significant-sites
summary (through :mod:`score_panels`) and writes, per test family present:

* ``score_heatmap_{family}.png`` -- rows = every MAG tested in any cell of the
  family; a phylum colour strip, then one tile per cell coloured by score
  (0-100 %, plasma) with the smallest q as stars; grey = untested (matplotlib);
* ``score_boxplot_{family}.png`` -- one box per cell, one dot per MAG, a line
  joining a MAG's two parallelism scores (plotnine);
* ``score_boxplot_{family}_by_phylum.png`` -- the same, one panel per phylum.

Plus ``score_figures.done`` listing the files written (the Snakemake sentinel).
Runs inside the pipeline after the scores, and standalone on any run directory::

    alleleflux-plot-scores --run_dir <run>/longitudinal \\
        --comparison pre_end-treatment_control --groups treatment control \\
        --output_dir <run>/longitudinal/figures/pre_end-treatment_control
"""

from __future__ import annotations

import argparse
import logging
import os

import matplotlib

matplotlib.use("Agg")  # headless: no display on a cluster node
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Patch
from plotnine import (
    aes,
    coord_cartesian,
    element_text,
    facet_wrap,
    geom_boxplot,
    geom_line,
    geom_point,
    geom_text,
    ggplot,
    labs,
    scale_x_continuous,
    theme,
    theme_bw,
)

from alleleflux.scripts.utilities.logging_config import setup_logging
from alleleflux.scripts.visualization.score_panels import (
    FAMILIES,
    FamilyData,
    load_family,
    q_stars,
    resolve_family,
    taxonomy_order,
    tree_order,
)

logger = logging.getLogger(__name__)

UNTESTED_GREY = "#d9d9d9"  # a cell the MAG was not tested in: never 0, never blank
SPECIES_WIDTH = 40  # row label = MAG id + species, species cut to this many characters
ROW_INCH = 0.15  # figure height grows with the row count so eighty rows stay legible


def _column_labels(data: FamilyData) -> list[str]:
    """Cell labels, flagging a cell the pipeline wrote no table for."""
    return [
        c.label if c.score_file else f"{c.label}\n(no eligible MAGs)"
        for c in data.family.cells
    ]


def draw_score_heatmap(
    data: FamilyData,
    order: list[str],
    absent_from_tree: list[str],
    title: str,
    out_path: str,
    q_threshold: float,
) -> Figure:
    """[row labels] [phylum strip] [score tiles with q stars]  + colour bar + phylum legend.

    Tiles: score 0-100 on plasma, grey where the MAG was not tested, the smallest
    q as white stars.  Rows: ``order`` then ``absent_from_tree`` (tree order
    only), the two blocks separated by a thin line.  The strip is a one-column
    image, so it scales with the rows like the tiles do.
    """
    rows = list(order) + list(absent_from_tree)
    scores = data.scores.reindex(rows)
    min_q = data.min_q.reindex(rows)
    tax = data.taxonomy.reindex(rows)
    n_rows, n_cols = len(rows), scores.shape[1]

    fig, (ax_phy, ax) = plt.subplots(
        1,
        2,
        figsize=(
            8.5,
            max(3.0, ROW_INCH * n_rows + 1.5),
        ),  # labels left, bar + legend right
        width_ratios=[0.3, n_cols],
        layout="constrained",  # keeps the colour bar, legend and title from colliding
    )
    # phylum strip: sorted phyla onto tab20 (cycling past 20), unknown phylum grey
    phylum = tax["phylum"].where(tax["phylum"].apply(lambda v: isinstance(v, str)))
    phyla = sorted(phylum.dropna().unique())
    colours = [matplotlib.colormaps["tab20"](i % 20) for i in range(len(phyla))]
    strip = np.array(
        [[phyla.index(v)] if isinstance(v, str) else [np.nan] for v in phylum]
    )
    strip_cmap = ListedColormap(colours or [UNTESTED_GREY])
    strip_cmap.set_bad(UNTESTED_GREY)
    ax_phy.imshow(
        np.ma.masked_invalid(strip),
        cmap=strip_cmap,
        aspect="auto",
        vmin=-0.5,
        vmax=max(len(phyla), 1) - 0.5,  # one colour bin per phylum
        interpolation="nearest",
    )
    ax_phy.set_xticks([])
    ax_phy.set_title("phylum", fontsize=8)
    # row labels (MAG id + species) hang off the strip, i.e. at the far left
    species = tax["species"].fillna("").astype(str).str.slice(0, SPECIES_WIDTH)
    ax_phy.set_yticks(range(n_rows))
    ax_phy.set_yticklabels([f"{m}  {sp}" for m, sp in zip(rows, species)], fontsize=6)
    ax_phy.tick_params(axis="y", length=0)

    # score tiles: grey where untested, never 0
    cmap = matplotlib.colormaps["plasma"].copy()
    cmap.set_bad(UNTESTED_GREY)
    im = ax.imshow(
        np.ma.masked_invalid(scores.to_numpy(dtype=float)),
        cmap=cmap,
        aspect="auto",
        vmin=0,
        vmax=100,
        interpolation="nearest",
    )
    for i in range(n_rows):  # q stars on the tiles; untested tiles have NaN q -> ""
        for j in range(n_cols):
            star = q_stars(min_q.iat[i, j], q_threshold)
            if star:
                ax.text(j, i, star, ha="center", va="center", fontsize=7, color="white")
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(_column_labels(data), rotation=45, ha="right", fontsize=7)
    ax.set_yticks([])
    if absent_from_tree:  # divider above the rows the tree did not contain
        for a in (ax_phy, ax):
            a.axhline(len(order) - 0.5, color="black", linewidth=0.8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.03, shrink=0.3, anchor=(0, 1))
    cbar.set_label("score (% sites significant)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    handles = [Patch(color=c, label=p) for c, p in zip(colours, phyla)]
    if phylum.isna().any():
        handles.append(Patch(color=UNTESTED_GREY, label="unknown"))
    ax.legend(
        handles=handles,
        title="phylum",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.62),  # below the colour bar, outside the tiles
        fontsize=6,
        title_fontsize=7,
        frameon=False,
    )
    fig.suptitle(title, fontsize=9)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")  # 6-pt labels need the pixels
    plt.close(fig)
    return fig


def _long_scores(data: FamilyData, rng) -> pd.DataFrame:
    """The score grid as one row per (MAG, cell) that was tested, ready for plotnine.

    Columns: ``mag``, ``cell`` (label), ``cell_idx`` (0, 1, 2 in family order),
    ``x`` (``cell_idx`` plus ONE horizontal jitter per MAG, reused for all its
    dots so a line between two of them lands exactly on both), ``score``,
    ``phylum``, ``paired`` (True on the two parallelism rows of a MAG that has
    both parallelism scores: those rows get the joining line).

    Example: M4 tested only in "treatment parallelism" -> one row, paired False.
    """
    cols = list(data.scores.columns)
    jitter = pd.Series(
        rng.uniform(-0.15, 0.15, len(data.scores)), index=data.scores.index
    )
    long = (
        data.scores.rename_axis("mag")
        .reset_index()
        .melt(id_vars="mag", var_name="cell", value_name="score")
        .dropna(subset=["score"])
    )
    long["cell_idx"] = long["cell"].map({c: i for i, c in enumerate(cols)})
    long["x"] = long["cell_idx"] + long["mag"].map(jitter)
    long["phylum"] = long["mag"].map(data.taxonomy["phylum"])
    par = [c.label for c in data.family.cells if c.kind == "parallelism"]
    both = data.scores[par].dropna().index if len(par) == 2 else []
    long["paired"] = long["mag"].isin(both) & long["cell"].isin(par)
    return long.reset_index(drop=True)


def _box_plot(long: pd.DataFrame, data: FamilyData, title: str, facet: bool):
    """One box per cell, a dot per MAG, a grey line between a MAG's two parallelism dots.

    ``n=`` per cell is written under each box (per panel when faceted).
    """
    cols = list(data.scores.columns)
    labels = _column_labels(data)
    group_cols = ["panel", "cell_idx"] if facet else ["cell_idx"]
    counts = long.groupby(group_cols, observed=True).size().reset_index(name="n")
    counts["label"] = "n=" + counts["n"].astype(str)
    # a cell with no tested MAG still gets its slot and an n=0
    frame = pd.DataFrame({"cell_idx": range(len(cols))})
    if facet:
        frame = frame.merge(
            pd.DataFrame({"panel": long["panel"].unique()}), how="cross"
        )
    counts = frame.merge(counts, how="left").fillna({"n": 0, "label": "n=0"})
    p = (
        ggplot(long)
        + geom_boxplot(
            aes(x="cell_idx", y="score", group="cell_idx"),
            width=0.5,
            outlier_alpha=0,  # every MAG is drawn as a dot anyway
            fill="white",
        )
        + geom_line(  # drawn under the dots
            aes(x="x", y="score", group="mag"),
            data=long[long["paired"]],
            color="grey",
            alpha=0.4,
            size=0.5,
        )
        + geom_point(aes(x="x", y="score", color="cell"), size=1.5, alpha=0.7)
        + geom_text(aes(x="cell_idx", y=-5, label="label"), data=counts, size=7)
        + scale_x_continuous(breaks=list(range(len(cols))), labels=labels)
        + coord_cartesian(xlim=(-0.6, len(cols) - 0.4), ylim=(-8, 102))
        + labs(x="", y="score (% sites significant)", title=title)
        + theme_bw()
        + theme(
            legend_position="none",
            axis_text_x=element_text(size=8),
            plot_title=element_text(size=10),
        )
    )
    if facet:
        panels = list(long["panel"].unique())
        p = (
            p
            + facet_wrap("~panel", nrow=1)
            + theme(
                figure_size=(len(panels) * (1.5 * len(cols) + 1), 4.5),
                strip_text=element_text(size=9),
            )
        )
    else:
        p = p + theme(figure_size=(1.8 * len(cols) + 2, 4.5))
    return p


def draw_score_boxplot(data: FamilyData, title: str, out_path: str):
    """One box per cell with every MAG as a dot; paired parallelism dots are joined."""
    long = _long_scores(data, np.random.default_rng(0))
    p = _box_plot(long, data, title, facet=False)
    p.save(out_path, dpi=150, verbose=False)
    return p


def draw_score_boxplot_by_phylum(
    data: FamilyData, title: str, out_path: str, min_mags_per_phylum: int
):
    """``draw_score_boxplot`` per phylum; phyla with fewer MAGs than the floor are pooled."""
    long = _long_scores(data, np.random.default_rng(0))
    phylum = data.taxonomy["phylum"].reindex(data.scores.index)
    counts = phylum.value_counts()
    big = sorted(counts[counts >= min_mags_per_phylum].index)
    small = sorted(set(counts.index) - set(big))
    # panel name per MAG: its phylum when the phylum has enough MAGs, else pooled
    n_pooled = int((~phylum.isin(big)).sum())
    panel_of = {
        m: (
            f"{p} (n={int(counts[p])})"
            if p in big
            else f"other ({len(small)} phyla, n={n_pooled})"
        )
        for m, p in phylum.items()
    }
    order = [f"{p} (n={int(counts[p])})" for p in big] + (
        [f"other ({len(small)} phyla, n={n_pooled})"] if n_pooled else []
    )
    long["panel"] = pd.Categorical(long["mag"].map(panel_of), categories=order)
    p = _box_plot(long, data, title, facet=True)
    p.save(out_path, dpi=150, verbose=False)
    return p


def paint_the_scoreboard(args: argparse.Namespace) -> list[str]:
    """Orchestrator: every family present -> three figures each; returns the files written."""
    scores_dir = os.path.join(args.run_dir, "scores", "processed", "combined", "MAG")
    cell_stats = pd.read_csv(
        os.path.join(
            args.run_dir,
            "p_value_summary",
            "significant_sites_summary",
            "significant_sites_mag_cell_stats_long.tsv",
        ),
        sep="\t",
        dtype={
            "mag_id": str,
            "group_analyzed": str,
            "test_type": str,
            "comparison": str,
        },
    )
    os.makedirs(args.output_dir, exist_ok=True)
    # PNG is always written: the run report embeds PNGs.  "pdf" adds a PDF copy
    # of every figure ("both" is the same thing, kept for older configs).
    formats = ["png"] + (["pdf"] if args.format in ("pdf", "both") else [])
    written: list[str] = []
    for family in args.families or FAMILIES:
        fam = resolve_family(
            scores_dir, args.comparison, args.groups, family, args.data_type
        )
        if fam is None:
            logger.info(
                f"{args.comparison}: no score table for family {family}: skipped"
            )
            continue
        data = load_family(fam, cell_stats, args.comparison)
        if args.data_type == "single":
            logger.info(f"{family}: single timepoint, divergence only")
        if args.tree:
            order, absent = tree_order(list(data.scores.index), args.tree)
            if absent:
                logger.info(
                    f"{family}: {len(absent)} MAG(s) absent from the tree drawn below a divider"
                )
        else:
            order, absent = taxonomy_order(data.taxonomy), []
        title = f"{args.comparison}  |  {family}" + (
            f"  |  {args.run_name}" if args.run_name else ""
        )
        for fmt in formats:
            for stem, draw in (
                (
                    f"score_heatmap_{family}",
                    lambda p: draw_score_heatmap(
                        data, order, absent, title, p, args.q_threshold
                    ),
                ),
                (
                    f"score_boxplot_{family}",
                    lambda p: draw_score_boxplot(data, title, p),
                ),
                (
                    f"score_boxplot_{family}_by_phylum",
                    lambda p: draw_score_boxplot_by_phylum(
                        data, title, p, args.min_mags_per_phylum
                    ),
                ),
            ):
                path = os.path.join(args.output_dir, f"{stem}.{fmt}")
                draw(path)
                written.append(path)
        logger.info(
            f"{family}: {len(data.scores)} MAGs, {data.scores.shape[1]} cells, figures written"
        )
    with open(os.path.join(args.output_dir, "score_figures.done"), "w") as fh:
        fh.write("\n".join(written) + ("\n" if written else ""))
    return written


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="Score heatmaps and boxplots per test family for one comparison of a finished run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--run_dir",
        required=True,
        help="the run's data-type directory (holds scores/ and p_value_summary/)",
    )
    parser.add_argument(
        "--comparison", required=True, help="e.g. pre_end-treatment_control"
    )
    parser.add_argument(
        "--groups",
        nargs=2,
        required=True,
        metavar=("GROUP_A", "GROUP_B"),
        help="the comparison's two group labels, in configuration order",
    )
    parser.add_argument("--output_dir", required=True)
    parser.add_argument(
        "--families",
        nargs="*",
        choices=FAMILIES,
        default=None,
        help=(
            "families to draw, e.g. '--families ttest' or '--families ttest lmm cmh'; "
            "ttest/wilcoxon use the paired divergence test, ttest_unpaired/mannwhitney "
            "the unpaired one (the one to read when samples cannot be paired across "
            "timepoints); default: every family with a score table"
        ),
    )
    parser.add_argument(
        "--tree",
        default=None,
        help="optional Newick tree with MAG ids as tips: rows follow its tip order",
    )
    parser.add_argument(
        "--format",
        choices=["png", "pdf", "both"],
        default="png",
        help="png: PNG only; pdf (or both): PNG plus a PDF copy of every figure",
    )
    parser.add_argument(
        "--q_threshold",
        type=float,
        default=0.05,
        help="q below which a cell gets a star (the run's p_value_threshold)",
    )
    parser.add_argument("--min_mags_per_phylum", type=int, default=3)
    parser.add_argument(
        "--data_type", choices=["longitudinal", "single"], default="longitudinal"
    )
    parser.add_argument("--run_name", default="", help="shown in figure titles")
    paint_the_scoreboard(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
