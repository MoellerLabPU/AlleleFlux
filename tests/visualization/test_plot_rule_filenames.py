"""The plotting rule's expected filenames must be the ones the plot script writes.

``plotting.smk`` (visualization workflow) predicts the output of
``alleleflux-plot-trajectories`` from the config alone::

    {mag}_combined_n{n_label}_{x_col_filename}_line{tag}.{fmt}

with ``n_label`` = "ALL" or the integer, ``x_col_filename`` = "bin_midpoint"
when day binning is on (else the configured x_col) and one ``tag`` per bin
width (``_bin<W>d``, or "" unbinned).  If that arithmetic drifts from what
``plot_allele_trajectory.py`` really names its files, Snakemake declares the
rule failed ("missing output") even though the script succeeded.

Two layers, so a mismatch is localised:

(a) ``TestRuleFilenameArithmetic`` execs the top-level parameter block of
    ``plotting.smk`` (everything above its first ``def``) with a stub
    ``config`` -- the same slice-and-exec trick as
    ``tests/test_site_depth_filter_params.py`` -- and checks the derived
    ``n_label`` / ``x_col_filename`` / ``bin_tags`` for three configs.
(b) ``TestScriptWritesWhatTheRuleExpects`` runs the real plot script as a
    subprocess on a tiny long table shaped exactly like the tracker's
    ``<mag>_frequency_table.long.tsv`` and asserts the output directory holds
    precisely the files the rule's arithmetic predicts (no misses, no extras).

Fixture geometry for (b): 2 mice (one per group) x 2 timepoints (day 0 and
day 25) x 3 anchor sites = 12 rows.  Days 0 and 25 land in different bins for
both widths used (10 -> "0-9"/"20-29", 20 -> "0-19"/"20-39"), so every binned
line plot still has two x-values to connect.

The three configs:

    A  n_sites_line 2,     unbinned            -> MAG_A_combined_n2_day_line.png
    B  n_sites_line "all", bin_width_days 10   -> MAG_A_combined_nALL_bin_midpoint_line_bin10d.png
    C  n_sites_line 2,     bin_width_days [10, 20]
         -> MAG_A_combined_n2_bin_midpoint_line_bin10d.png
            MAG_A_combined_n2_bin_midpoint_line_bin20d.png
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
PLOTTING_SMK = (
    REPO / "alleleflux" / "smk_workflow" / "visualization" / "smks" / "plotting.smk"
)

MAG = "MAG_A"
X_COL = "day"
FMT = "png"
VALUE_COL = "min_p_value"

# The three configs under test, as the ``plotting_params`` block of the
# visualization config.  Only the keys the filename arithmetic reads are
# load-bearing: n_sites_line, x_col, output_format, bin_width_days.
CONFIGS = {
    "A_unbinned": {"n_sites_line": 2, "bin_width_days": None},
    "B_bin10": {"n_sites_line": "all", "bin_width_days": 10},
    "C_bin10_20": {"n_sites_line": 2, "bin_width_days": [10, 20]},
}

# What each config's rule arithmetic must produce.
EXPECTED_ARITHMETIC = {
    "A_unbinned": ("2", "day", [""]),
    "B_bin10": ("ALL", "bin_midpoint", ["_bin10d"]),
    "C_bin10_20": ("2", "bin_midpoint", ["_bin10d", "_bin20d"]),
}


def _plotting_params(n_sites_line, bin_width_days) -> dict:
    """A plotting_params block; ``bin_width_days`` None means 'key absent'."""
    block = {
        "value_col": VALUE_COL,
        "n_sites_line": n_sites_line,
        "x_col": X_COL,
        "x_order": [],
        "plot_types": ["line"],
        "per_site": False,
        "output_format": FMT,
    }
    if bin_width_days is not None:
        block["bin_width_days"] = bin_width_days
    return block


def _rule_arithmetic(plotting_params: dict) -> dict:
    """Exec the parameter block at the top of ``plotting.smk`` with a stub config.

    The block is everything above the first top-level ``def`` -- it reads only
    ``config`` and defines n_label, x_col, output_fmt, x_col_filename,
    bin_widths, bin_tags and combined_per_replicate.  ``expand`` / checkpoints
    live below the slice and are not needed for the filename arithmetic.
    """
    lines = PLOTTING_SMK.read_text().splitlines()
    first_def = next(i for i, ln in enumerate(lines) if ln.startswith("def "))
    src = "\n".join(lines[:first_def])
    namespace = {"config": {"plotting_params": plotting_params}}
    exec(compile(src, str(PLOTTING_SMK), "exec"), namespace)
    return namespace


def _rule_expected_filenames(plotting_params: dict, mag: str = MAG) -> list[str]:
    """The basenames the rule's ``sentinels`` output would expand to for one MAG.

    Mirrors the f-string in ``rule plot_trajectories`` / ``get_plotting_targets``:
    ``{mag}_combined_n{n_label}_{x_col_filename}_line{tag}.{output_fmt}``, one
    per entry of ``bin_tags``.
    """
    ns = _rule_arithmetic(plotting_params)
    return [
        f"{mag}_combined_n{ns['n_label']}_{ns['x_col_filename']}_line{tag}.{ns['output_fmt']}"
        for tag in ns["bin_tags"]
    ]


class TestRuleFilenameArithmetic:
    """(a) the rule's config -> name arithmetic, exec'd from the .smk source."""

    @pytest.mark.parametrize("name", list(CONFIGS))
    def test_derived_parameters(self, name):
        ns = _rule_arithmetic(_plotting_params(**CONFIGS[name]))
        n_label, x_col_filename, bin_tags = EXPECTED_ARITHMETIC[name]
        assert ns["n_label"] == n_label
        assert ns["x_col_filename"] == x_col_filename
        assert ns["bin_tags"] == bin_tags
        assert ns["output_fmt"] == FMT

    def test_bin_widths_list_vs_scalar(self):
        # a scalar width is promoted to a one-element list (back-compat path)
        assert _rule_arithmetic(_plotting_params(2, 10))["bin_widths"] == [10]
        assert _rule_arithmetic(_plotting_params(2, [10, 20]))["bin_widths"] == [
            10,
            20,
        ]
        assert _rule_arithmetic(_plotting_params(2, None))["bin_widths"] == []

    def test_expected_filenames(self):
        assert _rule_expected_filenames(_plotting_params(2, None)) == [
            f"{MAG}_combined_n2_day_line.png"
        ]
        assert _rule_expected_filenames(_plotting_params("all", 10)) == [
            f"{MAG}_combined_nALL_bin_midpoint_line_bin10d.png"
        ]
        assert _rule_expected_filenames(_plotting_params(2, [10, 20])) == [
            f"{MAG}_combined_n2_bin_midpoint_line_bin10d.png",
            f"{MAG}_combined_n2_bin_midpoint_line_bin20d.png",
        ]


# --------------------------------------------------------------------------
# (b) the real script, on a tracker-shaped long table
# --------------------------------------------------------------------------


def _cli() -> list[str]:
    """The console script from the SAME env as pytest, else ``python -m`` fallback."""
    exe = Path(sys.executable).parent / "alleleflux-plot-trajectories"
    if exe.exists():
        return [str(exe)]
    return [
        sys.executable,
        "-m",
        "alleleflux.scripts.visualization.plot_allele_trajectory",
    ]


def _tracker_shaped_long_table(path: Path) -> None:
    """Write a long table with exactly the columns ``assemble_and_save_results`` emits.

    Column order follows ``save_long_version``: id columns (``min_p_value``
    present because the anchor file had it), the melted ``sample_id`` /
    ``frequency`` pair, then the metadata merge (``group``, ``time``,
    ``subjectID``, and the optional ``day`` + ``replicate``).

    Geometry: 2 mice x 2 timepoints x 3 sites.  Frequencies rise from the
    first to the second timepoint for every site so the line plots have a
    non-flat trajectory; the exact numbers do not matter to the filename test.
    """
    sites = [  # (contig, position, anchor_allele, min_p_value)
        ("c1", 100, "A", 0.001),
        ("c1", 200, "G", 0.010),
        ("c2", 50, "T", 0.020),
    ]
    samples = [  # (sample_id, group, time, subjectID, day, replicate)
        ("s1", "gA", "pre", "m1", 0, "r1"),
        ("s2", "gB", "pre", "m2", 0, "r2"),
        ("s3", "gA", "post", "m1", 25, "r1"),
        ("s4", "gB", "post", "m2", 25, "r2"),
    ]
    rows = []
    for contig, pos, allele, pval in sites:
        for k, (sid, group, time, subj, day, rep) in enumerate(samples):
            rows.append(
                {
                    "mag_id": MAG,
                    "contig": contig,
                    "position": pos,
                    "gene_id": f"gene_{contig}",
                    "anchor_allele": allele,
                    "anchor_criteria": "terminal_nucleotide_mean_freq",
                    "min_p_value": pval,
                    "sample_id": sid,
                    # 0.1/0.2 at pre, 0.6/0.7 at post: distinct per mouse and rising
                    "frequency": 0.1 * (k % 2 + 1) + (0.5 if day > 0 else 0.0),
                    "group": group,
                    "time": time,
                    "subjectID": subj,
                    "day": day,
                    "replicate": rep,
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def _run_plot_script(long_table: Path, out_dir: Path, plotting_params: dict):
    """Invoke the script the way ``rule plot_trajectories``'s shell block does.

    Only the flags the three configs need are passed (value_col, n_sites_line,
    x_col, plot_types, output_format, bin_width_days); the rule emits the
    other flags empty when their config keys are unset.
    """
    cmd = _cli() + [
        "--input_file", str(long_table),
        "--output_dir", str(out_dir),
        "--value_col", plotting_params["value_col"],
        "--n_sites_line", str(plotting_params["n_sites_line"]),
        "--x_col", plotting_params["x_col"],
        "--plot_types", *plotting_params["plot_types"],
        "--output_format", plotting_params["output_format"],
    ]  # fmt: skip
    widths = plotting_params.get("bin_width_days")
    if widths:
        widths = widths if isinstance(widths, list) else [widths]
        # the rule joins the widths into one --bin_width_days W [W ...]
        cmd += ["--bin_width_days", *map(str, widths)]
    env = {**os.environ, "MPLBACKEND": "Agg"}  # headless matplotlib
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


@pytest.fixture(scope="module")
def long_table(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("plot_rule") / f"{MAG}_frequency_table.long.tsv"
    _tracker_shaped_long_table(path)
    return path


class TestScriptWritesWhatTheRuleExpects:
    """(b) the script's real output file names, config by config."""

    @pytest.mark.parametrize("name", list(CONFIGS))
    def test_output_dir_holds_exactly_the_rule_sentinels(
        self, name, long_table, tmp_path
    ):
        params = _plotting_params(**CONFIGS[name])
        out_dir = tmp_path / "plots" / MAG  # the rule's params.out_dir layout
        done = _run_plot_script(long_table, out_dir, params)
        assert done.returncode == 0, f"stderr:\n{done.stderr}"

        expected = sorted(_rule_expected_filenames(params))
        actual = sorted(p.name for p in out_dir.iterdir())
        # Exact set equality: a renamed file would show up as both a missing
        # expectation and an unexpected extra, which pins the mismatch.
        assert actual == expected, (
            f"[{name}] rule expects {expected} but the script wrote {actual}\n"
            f"stderr:\n{done.stderr}"
        )
        for fname in expected:
            assert (out_dir / fname).stat().st_size > 0, f"{fname} is empty"
