"""CMH under a permuted (null) run on single-timepoint data.

A null run reuses the real run's tables and relabels ``group`` in memory from
the permuted sheet.  On single-timepoint data the CMH input is the
preprocessed between-groups table, which carries ``subjectID`` and
``replicate`` but NO ``sample_id`` -- the per-sample key the relabel helper
joins on by default.  The CLI must notice and join on ``subjectID`` instead
(one sample per subject on single data), not demand a column that is not
there.  Caught by the regression's permutation-equivalence test on 2026-10-30
(``ValueError: Usecols do not match columns ... ['sample_id']``).
"""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import pandas as pd

NUC = ["A", "C", "G", "T"]


def _single_preprocessed_table(path: Path) -> None:
    """Eight samples (4 control, 4 treatment, replicates A-D shared), 6 positions,
    with treatment carrying a T minority at every position so CMH has signal."""
    rows = []
    for pos in range(6):
        for subj, group, rep in [
            (1, "control", "A"),
            (2, "control", "B"),
            (3, "control", "C"),
            (4, "control", "D"),
            (5, "treatment", "A"),
            (6, "treatment", "B"),
            (7, "treatment", "C"),
            (8, "treatment", "D"),
        ]:
            a = 40 if group == "control" else 25
            t = 0 if group == "control" else 15
            cov = a + t
            rows.append(
                {
                    "contig": "c1",
                    "position": pos,
                    "gene_id": "g1",
                    "total_coverage": cov,
                    "A": a,
                    "C": 0,
                    "G": 0,
                    "T": t,
                    "group": group,
                    "subjectID": f"subj{subj}",
                    "replicate": rep,
                    **{
                        f"{n}_frequency": {"A": a, "C": 0, "G": 0, "T": t}[n] / cov
                        for n in NUC
                    },
                }
            )
    with gzip.open(path, "wt") as fh:
        pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)


def _permuted_sheet(path: Path) -> None:
    """The run's sample sheet with two subjects' groups swapped (subj3 <-> subj7)."""
    group_of = {
        1: "control",
        2: "control",
        3: "treatment",
        4: "control",
        5: "treatment",
        6: "treatment",
        7: "control",
        8: "treatment",
    }
    rep_of = {1: "A", 2: "B", 3: "C", 4: "D", 5: "A", 6: "B", 7: "C", 8: "D"}
    lines = ["sample_id\tbam_path\tsubjectID\tgroup\treplicate\ttime"]
    for s in range(1, 9):
        lines.append(
            f"sample{s}\tbams/sample{s}.bam\tsubj{s}\t{group_of[s]}\t{rep_of[s]}\tt1"
        )
    path.write_text("\n".join(lines) + "\n")


def _run_cmh(input_df: Path, out_dir: Path, permuted: Path | None):
    cmd = [
        "alleleflux-cmh",
        "--input_df",
        str(input_df),
        "--min_sample_num",
        "2",
        "--mag_id",
        "MAG_T",
        "--data_type",
        "single",
        "--groups",
        "treatment",
        "control",
        "--cpus",
        "1",
        "--output_dir",
        str(out_dir),
    ]
    if permuted is not None:
        cmd += ["--permuted_metadata", str(permuted)]
    return subprocess.run(cmd, capture_output=True, text=True)


def test_single_data_without_sample_id_relabels_by_subject(tmp_path):
    table = tmp_path / "MAG_T_allele_frequency_preprocessed.tsv.gz"
    _single_preprocessed_table(table)
    sheet = tmp_path / "permuted_metadata_treatment_control.tsv"
    _permuted_sheet(sheet)

    plain = _run_cmh(table, tmp_path / "plain", None)
    assert plain.returncode == 0, plain.stderr[-2000:]
    permuted = _run_cmh(table, tmp_path / "permuted", sheet)
    assert permuted.returncode == 0, permuted.stderr[-2000:]

    out_plain = pd.read_csv(tmp_path / "plain" / "MAG_T_cmh.tsv.gz", sep="\t")
    out_perm = pd.read_csv(tmp_path / "permuted" / "MAG_T_cmh.tsv.gz", sep="\t")
    assert len(out_plain) == 6 and len(out_perm) == 6
    pcol = next(c for c in out_perm.columns if "p_value" in c)
    # Swapping two subjects' labels must change the stratified test's answer:
    # the relabel really happened, by subjectID.
    assert not out_plain[pcol].round(8).equals(out_perm[pcol].round(8)), (
        out_plain[pcol].tolist(),
        out_perm[pcol].tolist(),
    )
