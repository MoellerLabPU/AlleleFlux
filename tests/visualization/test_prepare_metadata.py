"""alleleflux-prepare-metadata: rename a run's columns to the standard names and
combine several runs into one table, with no per-sample profile path (the tools
downstream read the allele-frequency cache)."""

import subprocess
import sys
from pathlib import Path

import pandas as pd


def _cli():
    exe = Path(sys.executable).parent / "alleleflux-prepare-metadata"
    return (
        [str(exe)]
        if exe.exists()
        else [
            sys.executable,
            "-m",
            "alleleflux.scripts.visualization.accessory.prepare_metadata",
        ]
    )


def _run(metadata_in, metadata_out, **cols):
    cmd = _cli() + [
        "--metadata-in",
        str(metadata_in),
        "--metadata-out",
        str(metadata_out),
    ]
    for flag, value in cols.items():
        cmd += [f"--{flag.replace('_', '-')}", value]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_renames_and_combines_two_runs_without_a_profile_column(tmp_path):
    # run A names its columns differently; run B uses the standard names
    a = tmp_path / "a.tsv"
    pd.DataFrame(
        {
            "sample": ["s1", "s2"],
            "diet": ["x", "y"],
            "tp": ["pre", "pre"],
            "mouse": ["m1", "m2"],
        }
    ).to_csv(a, sep="\t", index=False)
    b = tmp_path / "b.tsv"
    pd.DataFrame(
        {"sample_id": ["s3"], "group": ["x"], "time": ["post"], "subjectID": ["m1"]}
    ).to_csv(b, sep="\t", index=False)
    out = tmp_path / "combined.tsv"

    done = _run(
        a,
        out,
        sample_col="sample",
        group_col="diet",
        time_col="tp",
        subject_col="mouse",
    )
    assert done.returncode == 0, done.stderr
    done = _run(b, out)  # appends
    assert done.returncode == 0, done.stderr

    got = pd.read_csv(out, sep="\t")
    assert list(got.columns) == ["sample_id", "group", "time", "subjectID"]
    assert got["sample_id"].tolist() == ["s1", "s2", "s3"]
    assert got.loc[got.sample_id == "s1", "group"].item() == "x"
    assert "sample_profile_dir" not in got.columns


def test_duplicate_sample_across_runs_is_refused(tmp_path):
    a = tmp_path / "a.tsv"
    pd.DataFrame(
        {"sample_id": ["s1"], "group": ["x"], "time": ["pre"], "subjectID": ["m1"]}
    ).to_csv(a, sep="\t", index=False)
    out = tmp_path / "combined.tsv"
    assert _run(a, out).returncode == 0
    done = _run(a, out)
    assert done.returncode != 0 and "duplicate" in (done.stderr + done.stdout).lower()
