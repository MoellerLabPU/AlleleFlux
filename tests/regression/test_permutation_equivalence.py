"""Equivalence regression test for reuse-based permuted (null) runs.

The entire reuse-based permutation feature rests on one equality:

    reuse + relabel  ==  full rebuild from permuted metadata

This test proves it end-to-end on the bundled longitudinal example, computing
**both sides live** (no committed golden to maintain):

  1. **normal**  — a real run; its profiles/QC/allele-freq cache become the
     reuse source.
  2. **Path A (ground truth)** — the *old way*: a full pipeline whose
     ``metadata_path`` IS the permuted sheet, so the cache is rebuilt from the
     permuted labels.
  3. **Path B (new way)** — ``input.reuse_from`` the normal run +
     ``permutation`` enabled, so the cache is reused and ``group`` is relabeled
     in memory.

The statistical outputs of B must equal A within the same float32-aware scope
the e2e uses (``_scope.STAT_*``).  A sanity assertion also confirms the
permutation actually *changed* the result vs the normal run (so a no-op
permutation can't make the gate pass trivially).

Marked ``regression`` (≈3 pipeline runs, ~2-3 min); deselected from the fast
unit loop via ``-m "not regression"``.  Requires the ``alleleflux`` CLI in the
active conda env plus ``bwa``/``samtools`` (same prerequisites as the e2e).
"""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pytest
import yaml

from tests.regression.compare import compare_trees
from tests.regression.conftest import SCENARIOS
from tests.regression._scope import (
    STAT_ALLOW_GOLDEN_ONLY,
    STAT_ATOL,
    STAT_COMPARE_ONLY,
    STAT_RTOL,
)

logger = logging.getLogger(__name__)


def _run(cmd, cwd):
    """Run a pipeline command, failing the test with tails on non-zero exit."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(
            f"command failed ({' '.join(cmd)}):\n"
            f"  exit={result.returncode}\n"
            f"  stdout (tail):\n{result.stdout[-2000:]}\n"
            f"  stderr (tail):\n{result.stderr[-2000:]}"
        )
    return result


@pytest.mark.regression
def test_permutation_reuse_equals_rebuild(tmp_path):
    """reuse+relabel produces identical statistical outputs to a full rebuild."""
    bundle = SCENARIOS["longitudinal"]
    example_dir = Path(bundle["example_data_dir"])
    base_config = example_dir / bundle["config_name"]
    assert base_config.exists(), f"base config missing: {base_config}"

    with base_config.open() as f:
        cfg = yaml.safe_load(f)

    # All three runs use cwd=example_dir so the metadata's relative bam_path
    # ("bams/…") resolves; outputs go to absolute tmp_path dirs to avoid
    # polluting the committed example tree.
    normal_out = tmp_path / "normal"
    rebuild_out = tmp_path / "rebuild"
    reuse_out = tmp_path / "reuse"
    perm_dir = tmp_path / "permuted_metadata"
    perm_dir.mkdir(parents=True, exist_ok=True)

    def write_cfg(name: str, mutate) -> Path:
        c = yaml.safe_load(yaml.safe_dump(cfg))  # deep copy
        mutate(c)
        p = tmp_path / name
        with p.open("w") as fh:
            yaml.safe_dump(c, fh, sort_keys=False)
        return p

    # 1. Normal run — the reuse source.
    normal_cfg = write_cfg(
        "config_normal.yml",
        lambda c: c["output"].__setitem__("root_dir", str(normal_out)),
    )
    _run(["alleleflux", "run", "--config", str(normal_cfg), "--threads", "4"], example_dir)
    normal_data = normal_out / bundle["pipeline_data_dir"]
    assert normal_data.exists(), f"normal run produced no output at {normal_data}"

    # 2. Generate ONE permuted sheet for the treatment/control pair.  The example
    #    is a mixed-group-replicate design, so --block replicate keeps each
    #    replicate balanced (eligibility preserved → a non-trivial, non-empty
    #    comparison).  Named for the groups label the rules expect.
    treatment = cfg["analysis"]["groups_combinations"][0]["treatment"]
    control = cfg["analysis"]["groups_combinations"][0]["control"]
    perm_sheet = perm_dir / f"permuted_metadata_{treatment}_{control}.tsv"
    _run(
        [
            "alleleflux-permute-metadata",
            "--input", cfg["input"]["metadata_path"],
            "--output", str(perm_sheet),
            "--unit", "subjectID",
            "--block", "replicate",
            "--swaps", "1",
            "--groups", treatment, control,
            "--seed", "5",
        ],
        example_dir,
    )
    assert perm_sheet.exists()

    # 3. Path A — full rebuild from the permuted metadata (ground truth).
    rebuild_cfg = write_cfg(
        "config_rebuild.yml",
        lambda c: (
            c["input"].__setitem__("metadata_path", str(perm_sheet)),
            c["output"].__setitem__("root_dir", str(rebuild_out)),
        ),
    )
    _run(["alleleflux", "run", "--config", str(rebuild_cfg), "--threads", "4"], example_dir)
    rebuild_data = rebuild_out / bundle["pipeline_data_dir"]

    # 4. Path B — reuse normal cache + relabel.  A *leaf* config (permuted_metadata_dir
    #    set) so execute_workflow runs the pipeline directly rather than re-orchestrating.
    reuse_cfg = write_cfg(
        "config_reuse.yml",
        lambda c: (
            c["input"].__setitem__("reuse_from", str(normal_data)),
            c["output"].__setitem__("root_dir", str(reuse_out)),
            c.__setitem__(
                "permutation",
                {"enabled": True, "permuted_metadata_dir": str(perm_dir)},
            ),
        ),
    )
    _run(["alleleflux", "run", "--config", str(reuse_cfg), "--threads", "4"], example_dir)
    reuse_data = reuse_out / bundle["pipeline_data_dir"]

    # GATE: reuse (B) must equal rebuild (A) on the scoped statistical outputs.
    gate = compare_trees(
        reuse_data, rebuild_data,
        include_prefixes=STAT_COMPARE_ONLY,
        allow_golden_only=STAT_ALLOW_GOLDEN_ONLY,
        rtol=STAT_RTOL, atol=STAT_ATOL,
    )
    assert not gate, (
        f"reuse+relabel diverged from full rebuild ({len(gate)} file(s)):\n\n"
        + "\n\n".join(gate)
    )

    # SANITY: the permutation must actually change the result vs the normal run,
    # otherwise the gate could pass on a no-op permutation.
    sanity = compare_trees(
        rebuild_data, normal_data,
        include_prefixes=STAT_COMPARE_ONLY,
        allow_golden_only=STAT_ALLOW_GOLDEN_ONLY,
        rtol=STAT_RTOL, atol=STAT_ATOL,
    )
    assert sanity, (
        "permutation produced identical stats to the normal run — the gate "
        "would pass trivially; pick a seed/unit that actually relabels."
    )
