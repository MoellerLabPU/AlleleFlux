"""End-to-end test of the visualization workflow (``smk_workflow/visualization``).

Runs the real Snakemake workflow -- metadata prep -> terminal nucleotide
analysis (checkpoint) -> allele tracking -> plotting -- on a tiny synthetic
allele-frequency cache, twice: once unbinned on the ``time`` axis and once
binned by ``day``, so the plotting rule's expected filenames are checked
against what the plot script really writes.  Every table the workflow
produces is opened and its planted numbers asserted, not just its existence.

Why this exists: on 2026-10-19 a config key the metadata-prep rule still read
had been renamed, which would have crashed the first stage; no test ran the
workflow, so nothing noticed.  This one would have.

Fixture geometry (one MAG ``MAG_X``, one contig ``c1``, 2 mice x 2 timepoints,
all in group ``g1``; counts are A,C,G,T):

    site  m1_pre         m2_pre         m1_post        m2_post        expected
    0     20,0,0,0       20,0,0,0       2,0,18,0       5,0,15,0       G by every method (rose 0 -> 0.9 / 0.75)
    1     0,0,20,0       0,0,20,0       0,0,20,0       0,0,20,0       G by mean/majority; freq_change blank (no change)
    2     10,0,10,0      12,0,8,0       16,0,4,0       18,0,2,0       A by every method (rose 0.5/0.6 -> 0.8/0.9)

About 30 s per configuration (the plot script dominates).  Lives under
tests/regression/ and carries the ``regression`` marker, so both the folder
habit (``--ignore=tests/regression``) and the label habit (``-m "not regression"``)
leave it out of the fast unit loop; ``pytest tests/regression/`` runs it with the
pipeline golden test.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from tests.utilities.test_cache_io import write_cache

pytestmark = pytest.mark.regression  # slow: runs a real Snakemake workflow (~50 s)

SNAKEFILE = (
    Path(__file__).resolve().parents[3]
    / "alleleflux"
    / "smk_workflow"
    / "visualization"
    / "visualization.smk"
)

PRE = {
    "m1": [(20, 0, 0, 0), (0, 0, 20, 0), (10, 0, 10, 0)],
    "m2": [(20, 0, 0, 0), (0, 0, 20, 0), (12, 0, 8, 0)],
}
POST = {
    "m1": [(2, 0, 18, 0), (0, 0, 20, 0), (16, 0, 4, 0)],
    "m2": [(5, 0, 15, 0), (0, 0, 20, 0), (18, 0, 2, 0)],
}


def _write_inputs(root: Path) -> dict:
    """The cache, the metadata and the significant-sites table under ``root``."""
    cache = root / "allele_freq_cache"
    for tp, counts in (("pre", PRE), ("post", POST)):
        rows = [
            ("c1", pos, f"{mouse}_{tp}", "g1", mouse, *cnt)
            for mouse, sites in counts.items()
            for pos, cnt in enumerate(sites)
        ]
        write_cache(str(cache), "MAG_X", tp, rows)
    metadata = root / "metadata.tsv"
    pd.DataFrame(
        [
            {
                "sample_id": f"{mouse}_{tp}",
                "subjectID": mouse,
                "group": "g1",
                "time": tp,
                "day": 0 if tp == "pre" else 30,  # for the binned configuration
            }
            for mouse in ("m1", "m2")
            for tp in ("pre", "post")
        ]
    ).to_csv(metadata, sep="\t", index=False)
    sig = root / "sig.tsv"
    pd.DataFrame(
        {
            "period": "pre_post",
            "mag_id": "MAG_X",
            "contig": "c1",
            "position": [0, 1, 2],
            "gene_id": ["", "", ""],
            "test_type": "two_sample_paired_tTest",
            "min_p_value": 0.0,
            "source_file": "x",
            "q_value": 0.0,
        }
    ).to_csv(sig, sep="\t", index=False)
    return {"cache": cache, "metadata": metadata, "sig": sig}


def _config(inputs: dict, out: Path, binned: bool) -> dict:
    """A complete visualization config, mirroring visualization_config.yaml's keys."""
    plotting = {
        "value_col": "min_p_value",
        "n_sites_line": 3,
        "n_sites_per_site": 3,
        "x_col": "day" if binned else "time",
        "x_order": [] if binned else ["pre", "post"],
        "plot_types": ["line"],
        "per_site": False,
        "output_format": "png",
    }
    if binned:
        plotting["bin_width_days"] = 30
        plotting["min_samples_per_bin"] = 1
    return {
        "run_name": "e2e",
        "output_dir": str(out),
        "resources": {"threads_per_job": 1, "mem_per_job": "1G"},
        "metadata_params": [{"name": "run_a", "path": str(inputs["metadata"])}],
        "significant_sites_file": str(inputs["sig"]),
        "cache_dir": str(inputs["cache"]),
        "metadata": str(inputs["metadata"]),
        "analysis_params": {
            "group": "g1",
            "timepoint": "post",
            "p_value_column": "q_value",
            "p_value_threshold": 0.05,
            "test_type": "two_sample_paired_tTest",
            "initial_timepoint": "pre",
        },
        "tracking_params": {
            "anchor_column": "terminal_nucleotide_majority_vote",
            "min_cov_per_site": 3,
        },
        "plotting_params": plotting,
    }


def _snakemake(configfile: Path, cwd: Path, *extra: str) -> subprocess.CompletedProcess:
    exe = Path(sys.executable).parent / "snakemake"
    cmd = [
        str(exe) if exe.exists() else "snakemake",
        "--snakefile",
        str(SNAKEFILE),
        "--configfile",
        str(configfile),
        "--cores",
        "2",
        *extra,
    ]
    env = {**os.environ, "MPLBACKEND": "Agg"}
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=env)


@pytest.mark.parametrize("binned", [False, True], ids=["time_axis", "binned_by_day"])
def test_workflow_runs_end_to_end_on_a_synthetic_cache(tmp_path, binned):
    inputs = _write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    out.mkdir()
    configfile = tmp_path / "config.yaml"
    configfile.write_text(yaml.safe_dump(_config(inputs, out, binned), sort_keys=False))

    # ---- dry run: the pre-checkpoint DAG resolves (metadata prep + terminal analysis)
    dry = _snakemake(configfile, out, "--dry-run")
    assert dry.returncode == 0, dry.stderr + dry.stdout
    for rule in ("clean_metadata", "aggregate_metadata", "terminal_analysis"):
        assert rule in dry.stdout, f"{rule} not scheduled:\n{dry.stdout}"

    # ---- real run
    done = _snakemake(configfile, out)
    assert done.returncode == 0, done.stderr + done.stdout
    assert "steps (100%) done" in (done.stdout + done.stderr)

    # stage 1: the combined metadata, standard columns only
    combined = pd.read_csv(out / "combined_metadata.tsv", sep="\t")
    assert sorted(combined["sample_id"]) == ["m1_post", "m1_pre", "m2_post", "m2_pre"]
    assert "sample_profile_dir" not in combined.columns

    # stage 2: the terminal nucleotides, per method
    term = pd.read_csv(
        out / "terminal_nucleotide" / "MAG_X" / "MAG_X_terminal_nucleotides.tsv",
        sep="\t",
    ).set_index("position")
    assert term["terminal_nucleotide_mean_freq"].tolist() == ["G", "G", "A"]
    assert term["terminal_nucleotide_majority_vote"].tolist() == ["G", "G", "A"]
    change = term["terminal_nucleotide_freq_change"]
    assert change[0] == "G" and change[2] == "A" and pd.isna(change[1])
    summary = pd.read_csv(
        out / "terminal_nucleotide" / "terminal_nucleotide_analysis_summary.tsv",
        sep="\t",
    )
    assert summary["mag_id"].tolist() == ["MAG_X"]
    assert int(summary["sites_processed"].iloc[0]) == 3

    # stage 3: the tracked frequencies of the majority-vote anchor, per sample
    wide = pd.read_csv(
        out / "track_freqs" / "MAG_X_frequency_table.tsv", sep="\t"
    ).set_index("position")
    np.testing.assert_allclose(
        wide.loc[0, ["m1_pre", "m2_pre", "m1_post", "m2_post"]].astype(float),
        [0.0, 0.0, 0.9, 0.75],
    )
    np.testing.assert_allclose(
        wide.loc[2, ["m1_pre", "m2_pre", "m1_post", "m2_post"]].astype(float),
        [0.5, 0.6, 0.8, 0.9],
    )
    long = pd.read_csv(out / "track_freqs" / "MAG_X_frequency_table.long.tsv", sep="\t")
    assert len(long) == 3 * 4  # sites x samples
    assert set(long["subjectID"]) == {"m1", "m2"}

    # stage 4: the plot the rule declares as its sentinel exists, with the
    # filename arithmetic of plotting.smk (unbinned: x_col; binned: bin_midpoint + tag)
    expected = (
        "MAG_X_combined_n3_bin_midpoint_line_bin30d.png"
        if binned
        else "MAG_X_combined_n3_time_line.png"
    )
    plot = out / "plots" / "MAG_X" / expected
    assert plot.exists(), sorted(p.name for p in (out / "plots" / "MAG_X").glob("*"))
    assert plot.stat().st_size > 0


def test_rerun_is_a_no_op(tmp_path):
    """A second invocation finds everything up to date (no rule reruns itself)."""
    inputs = _write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    out.mkdir()
    configfile = tmp_path / "config.yaml"
    configfile.write_text(yaml.safe_dump(_config(inputs, out, False), sort_keys=False))
    assert _snakemake(configfile, out).returncode == 0
    again = _snakemake(configfile, out)
    assert again.returncode == 0
    assert "Nothing to be done" in (again.stdout + again.stderr)
