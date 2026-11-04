#!/usr/bin/env python3
"""One self-contained HTML page per run: ``alleleflux-run-report``.

Reads a finished run's bookkeeping tables and figures and writes a single HTML
file (inline CSS, images embedded, no scripts) with, in order:

1. the run: name, directory, date, version, code commit, the saved runtime
   configuration (``quality_control`` shown expanded);
2. per comparison, the sample-and-MAG funnel: samples in the sheet, MAGs with
   profiles, MAGs eligible per test per the eligibility tables, MAGs still
   eligible after preprocessing;
3. per comparison, per test cell: MAGs tested, MAGs with at least one site
   with q below the threshold (count and percent of tested);
4. the site-depth filter: one row per timepoint and a histogram of the
   fraction removed per MAG;
5. every figure under the figures directory, indexed and captioned.

The per-MAG metadata lives under the *reuse* directory when the
run reused profiles (``input.reuse_from``), so those directories are taken
from the flags, else from the newest ``config_runtime_*.yml`` in the run root
(the parent of the data-type directory), else assumed under the run directory.
A section whose input is missing says so in one line; the rest renders.

    alleleflux-run-report --run_dir <run>/longitudinal \\
        --figures_dir <run>/longitudinal/figures \\
        --output <run>/longitudinal/figures/run_report.html \\
        --comparisons pre_end-treatment_control
"""

from __future__ import annotations

import argparse
import base64
import datetime as dt
import glob
import io
import logging
import os
import subprocess
from importlib import metadata

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from jinja2 import Environment, PackageLoader, select_autoescape

from alleleflux.scripts.utilities.logging_config import setup_logging

logger = logging.getLogger(__name__)

CAPTIONS = {
    "heatmap": (
        "Rows: every MAG tested in any cell of the family; the first column is the "
        "phylum. Tile colour: score, the % of a MAG's sites with p below the "
        "threshold (0-100). Stars: the smallest BH q at any site, *** q < 0.001, "
        "** q < 0.01, * q < threshold. Grey: the MAG was not tested in that cell."
    ),
    "boxplot": (
        "One box per cell, one dot per MAG; a grey line joins the same MAG's two "
        "parallelism scores. n = MAGs in the box."
    ),
    "boxplot_by_phylum": "The same boxplots, one panel per phylum; small phyla are pooled.",
}


# --------------------------------------------------------------------------- #
# Inputs
# --------------------------------------------------------------------------- #
def newest_runtime_config(run_root: str) -> dict | None:
    """The newest ``config_runtime_*.yml`` in ``run_root`` (names sort by timestamp), or None."""
    paths = sorted(glob.glob(os.path.join(run_root, "config_runtime_*.yml")))
    if not paths:
        return None
    with open(paths[-1]) as fh:
        return yaml.safe_load(fh) or {}


def sample_counts(metadata_path: str) -> pd.DataFrame:
    """Samples per (group, time) in the sheet; ``time`` is "" when the sheet has none."""
    meta = pd.read_csv(metadata_path, sep="\t", dtype=str)
    if "time" not in meta.columns:
        meta["time"] = ""
    out = (
        meta.groupby(["group", "time"], sort=True).size().reset_index(name="n_samples")
    )
    return out


def mags_with_profiles(metadata_dir: str, timepoints: str) -> int | None:
    """Distinct MAGs in ``inputMetadata_{timepoints}/summary_mag_index.tsv``, or None."""
    path = os.path.join(
        metadata_dir, f"inputMetadata_{timepoints}", "summary_mag_index.tsv"
    )
    if not os.path.exists(path):
        return None
    return int(
        pd.read_csv(path, sep="\t", dtype=str, usecols=["mag_id"])["mag_id"].nunique()
    )


def _true_counts(path: str) -> dict[str, int]:
    """``{column: n_true}`` for every yes/no column of the table.

    Flag columns are recognised by their values (only True/False), not their
    names: the preprocessing tables also carry ``positions_eligible_<x>``, an
    integer count that must not be reported as a flag.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    return {
        c: int((df[c] == "True").sum())
        for c in df.columns
        if set(df[c].dropna().unique()) <= {"True", "False"}
    }


def eligibility_counts(path: str) -> dict[str, int] | None:
    """``{column: n_true}`` for every ``*_eligible*`` column of the eligibility table, or None."""
    return _true_counts(path) if os.path.exists(path) else None


def preprocessing_counts(directory: str, comparison: str) -> dict[str, int] | None:
    """The same over the between-groups and within-groups preprocessing tables, or None."""
    out: dict[str, int] = {}
    for kind in ("between_groups", "within_groups"):
        path = os.path.join(
            directory, f"preprocessing_eligibility_{kind}_{comparison}.tsv"
        )
        if os.path.exists(path):
            out.update(_true_counts(path))
    return out or None


def cell_significance(
    cell_stats: pd.DataFrame, comparison: str, q_threshold: float
) -> pd.DataFrame:
    """Per test cell: MAGs tested, MAGs with >= 1 site with q < threshold, and % of tested.

    A MAG has a site below the threshold exactly when its smallest q is below it.
    Example: three MAGs with min q 0.001, 0.2, 0.049 at 0.05 -> tested 3,
    significant 2, 66.7 %.
    """
    cs = cell_stats[cell_stats["comparison"] == comparison].copy()
    cs["group_analyzed"] = cs["group_analyzed"].fillna("").astype(str)
    cs["sig"] = cs["min_q_value"].astype(float) < q_threshold
    g = cs.groupby(["test_type", "group_analyzed"], sort=True)
    out = pd.DataFrame(
        {
            "n_tested": g["mag_id"].nunique(),
            "n_significant": g["sig"].sum().astype(int),
        }
    ).reset_index()
    out["pct_significant"] = (100 * out["n_significant"] / out["n_tested"]).round(1)
    return out[
        [
            "test_type",
            "group_analyzed",
            "n_tested",
            "n_significant",
            "pct_significant",
        ]
    ]


def site_filter_section(directory: str) -> dict | None:
    """One summary row per timepoint from the per-(sample, MAG) filter table, plus a
    histogram (PNG bytes) of the per-MAG fraction removed.

    The raw table has one row per sample and MAG (thousands on a real run), so
    the page shows totals: pairs, MAGs, sites covered and removed, the overall
    percentage removed, and the median / max percentage over rows.

    Example: rows (covered 100, removed 5) and (covered 200, removed 15) at one
    timepoint -> sample_mag_pairs_checked 2, sites_covered 300, sites_removed 20,
    pct_removed 6.67, median_row_pct 6.25, max_row_pct 7.5.
    """
    summary = os.path.join(directory, "site_depth_filter_summary.tsv")
    by_mag = os.path.join(directory, "site_depth_filter_summary_by_mag.tsv")
    if not os.path.exists(summary):
        return None
    rows = pd.read_csv(summary, sep="\t", dtype={"mag_id": str, "sample_id": str})
    rows["removed"] = rows["n_removed_low"] + rows["n_removed_high"]
    rows["row_fraction"] = rows["removed"] / rows["n_sites_covered"]
    per_tp = (
        rows.groupby("timepoint")
        .agg(
            sample_mag_pairs_checked=("mag_id", "size"),
            mags=("mag_id", "nunique"),
            sites_covered=("n_sites_covered", "sum"),
            sites_removed=("removed", "sum"),
            median_row_fraction=("row_fraction", "median"),
            max_row_fraction=("row_fraction", "max"),
        )
        .reset_index()
    )
    # percentages, not fractions: the page prints one decimal, and 1.2% must not show as 0.0
    per_tp["pct_removed"] = 100 * per_tp["sites_removed"] / per_tp["sites_covered"]
    per_tp["median_row_pct"] = 100 * per_tp.pop("median_row_fraction")
    per_tp["max_row_pct"] = 100 * per_tp.pop("max_row_fraction")
    per_tp = per_tp[
        [
            "timepoint",
            "sample_mag_pairs_checked",
            "mags",
            "sites_covered",
            "sites_removed",
            "pct_removed",
            "median_row_pct",
            "max_row_pct",
        ]
    ]
    section = {"summary": per_tp, "histogram_png": b""}
    if os.path.exists(by_mag):
        df = pd.read_csv(by_mag, sep="\t")
        frac = [c for c in df.columns if "fraction" in c]
        if frac:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.hist(df[frac[0]].dropna(), bins=30)
            ax.set_xlabel(frac[0])
            ax.set_ylabel("MAG x timepoint rows")
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
            plt.close(fig)
            section["histogram_png"] = buf.getvalue()
    return section


def collect_figures(figures_dir: str) -> list[dict]:
    """Every ``score_*.png`` under ``figures_dir/<comparison>/``, embedded as base64."""
    out = []
    for path in sorted(glob.glob(os.path.join(figures_dir, "*", "score_*.png"))):
        stem = os.path.basename(path)[
            len("score_") : -len(".png")
        ]  # heatmap_ttest, boxplot_ttest_by_phylum
        kind_word, rest = stem.split("_", 1)
        family = rest.replace("_by_phylum", "")
        kind = kind_word + ("_by_phylum" if rest.endswith("_by_phylum") else "")
        with open(path, "rb") as fh:
            png_b64 = base64.b64encode(fh.read()).decode()
        out.append(
            {
                "comparison": os.path.basename(os.path.dirname(path)),
                "family": family,
                "kind": kind,
                "path": path,
                "caption": CAPTIONS.get(kind, ""),
                "png_b64": png_b64,
            }
        )
    return out


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
def render_report(ctx: dict) -> str:
    """The HTML for ``ctx`` (see ``build_context`` for its keys)."""
    env = Environment(
        loader=PackageLoader("alleleflux.scripts.visualization", "templates"),
        autoescape=select_autoescape(["html"]),
    )
    # integers get thousands separators (site counts run to billions), floats one decimal
    env.filters["commas"] = lambda v: (
        f"{v:,}" if isinstance(v, (int, np.integer)) else v
    )
    env.filters["table"] = lambda df: df.to_html(
        index=False,
        border=0,
        na_rep="",
        float_format=lambda v: f"{v:.1f}",
        formatters={c: "{:,}".format for c in df.select_dtypes("integer").columns},
    )
    return env.get_template("run_report.html.j2").render(**ctx)


def _funnel(args, cfg: dict, comparison: str) -> list[tuple[str, object]]:
    """The funnel rows of one comparison as (label, value); value None = not available."""
    tp = comparison.split("-", 1)[0]
    rows: list[tuple[str, object]] = []
    sheet = args.metadata_path or cfg.get("input", {}).get("metadata_path")
    if sheet and os.path.exists(sheet):
        for r in sample_counts(sheet).itertuples(index=False):
            label = f"samples in the metadata sheet ({r.group}" + (
                f", {r.time})" if r.time else ")"
            )
            rows.append((label, int(r.n_samples)))
    else:
        rows.append(("samples in the metadata sheet", None))
    rows.append(
        (f"MAGs with profiles ({tp})", mags_with_profiles(args.metadata_dir, tp))
    )
    elig = eligibility_counts(
        os.path.join(args.run_dir, f"eligibility_table_{comparison}.tsv")
    )
    if elig is None:
        rows.append(("MAGs with enough samples (eligibility table)", None))
    for col, n in (elig or {}).items():
        rows.append((f"MAGs with enough samples {_eligibility_label(col)}", n))
    pre = preprocessing_counts(
        os.path.join(args.run_dir, "preprocessing_eligibility"), comparison
    )
    if pre is None:
        rows.append(("MAGs with enough sites after preprocessing", None))
    for col, n in (pre or {}).items():
        rows.append(
            (f"MAGs with enough sites after preprocessing {_eligibility_label(col)}", n)
        )
    return rows


def _eligibility_label(column: str) -> str:
    """Plain-language test name for an eligibility flag column.

    The workflow gates each test on one flag: the unpaired flag feeds the
    unpaired two-sample tests and LMM, the paired flag the paired two-sample
    tests and CMH, and ``single_sample_eligible_<group>`` every within-group
    test of that group.  Example: ``single_sample_eligible_control`` ->
    "for the within-group tests in control (parallelism)".
    """
    if column in ("unpaired_test_eligible", "two_sample_unpaired_eligible"):
        return "for the unpaired between-group tests and LMM (divergence)"
    if column in ("paired_test_eligible", "two_sample_paired_eligible"):
        return "for the paired between-group tests and CMH (divergence)"
    if column.startswith("single_sample_eligible_"):
        group = column[len("single_sample_eligible_") :]
        return f"for the within-group tests in {group} (parallelism)"
    return f"({column})"


def build_context(args: argparse.Namespace) -> dict:
    """Everything the template needs, from the run directory and the flags."""
    run_root = os.path.dirname(os.path.abspath(args.run_dir.rstrip("/")))
    cfg = (
        newest_runtime_config(run_root)
        if args.config is None
        else (yaml.safe_load(open(args.config)) or {})
    )
    cfg = cfg or {}
    reuse = cfg.get("input", {}).get("reuse_from") or args.run_dir
    args.metadata_dir = args.metadata_dir or os.path.join(reuse, "inputMetadata")
    perm = cfg.get("permutation", {}) or {}
    permuted = bool(
        args.permuted or perm.get("enabled") or perm.get("permuted_metadata_dir")
    )

    cs_path = os.path.join(
        args.run_dir,
        "p_value_summary",
        "significant_sites_summary",
        "significant_sites_mag_cell_stats_long.tsv",
    )
    cell_stats = (
        pd.read_csv(cs_path, sep="\t", dtype={"group_analyzed": str})
        if os.path.exists(cs_path)
        else None
    )
    comparisons = []
    for cmp_ in args.comparisons:
        comparisons.append(
            {
                "name": cmp_,
                "funnel": _funnel(args, cfg, cmp_),
                "significance": (
                    cell_significance(cell_stats, cmp_, args.q_threshold)
                    if cell_stats is not None
                    else None
                ),
            }
        )
    try:
        version = metadata.version("alleleflux")
    except metadata.PackageNotFoundError:
        version = ""
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        commit = subprocess.run(
            ["git", "-C", here, "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
        ).stdout.strip()
    except Exception:  # not a checkout, no git: the field stays blank
        commit = ""
    site = site_filter_section(os.path.join(args.run_dir, "site_depth_filter"))
    if site is not None:
        site["histogram_b64"] = (
            base64.b64encode(site["histogram_png"]).decode()
            if site["histogram_png"]
            else ""
        )
    return {
        "run_name": cfg.get("run_name", ""),
        "run_dir": os.path.abspath(args.run_dir),
        "date": dt.date.today().isoformat(),
        "version": version,
        "commit": commit,
        "config_text": yaml.safe_dump(cfg, sort_keys=False) if cfg else "",
        "qc_text": (
            yaml.safe_dump(cfg.get("quality_control", {}), sort_keys=False)
            if cfg
            else ""
        ),
        "permuted": permuted,
        "data_type": args.data_type,
        "comparisons": comparisons,
        "site_filter": site,
        "figures": collect_figures(args.figures_dir),
        "q_threshold": args.q_threshold,
    }


def compile_the_dossier(args: argparse.Namespace) -> str:
    """Orchestrator: context -> HTML -> ``args.output``; returns the path."""
    html = render_report(build_context(args))
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        fh.write(html)
    logger.info(f"wrote {args.output} ({len(html) / 1e6:.1f} MB)")
    return args.output


def main():
    setup_logging()
    parser = argparse.ArgumentParser(
        description="One self-contained HTML report for a finished run: funnel, significant MAGs per test, figures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--run_dir", required=True, help="the run's data-type directory"
    )
    parser.add_argument(
        "--figures_dir",
        required=True,
        help="where alleleflux-plot-scores wrote figures/<comparison>/",
    )
    parser.add_argument("--output", required=True, help="the HTML file to write")
    parser.add_argument(
        "--comparisons", nargs="+", required=True, help="e.g. pre_end-treatment_control"
    )
    parser.add_argument(
        "--metadata_dir",
        default=None,
        help="inputMetadata/ directory (default: <reuse_from or run_dir>/inputMetadata)",
    )
    parser.add_argument(
        "--metadata_path",
        default=None,
        help="the sample sheet (default: input.metadata_path of the runtime config)",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="runtime config YAML (default: newest config_runtime_*.yml in the run root)",
    )
    parser.add_argument("--q_threshold", type=float, default=0.05)
    parser.add_argument(
        "--data_type", choices=["longitudinal", "single"], default="longitudinal"
    )
    parser.add_argument(
        "--permuted",
        action="store_true",
        help="mark the report as a permuted (null) run",
    )
    compile_the_dossier(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
