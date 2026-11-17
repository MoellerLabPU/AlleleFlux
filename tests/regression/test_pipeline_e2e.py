"""
End-to-end regression test for the AlleleFlux pipeline.

Runs the bundled example dataset through the full Snakemake workflow and
compares every output file to a committed golden snapshot.  Catches any
change in pipeline behaviour that affects downstream data — including
silent numerical drift that unit tests cannot see.

Siblings in this folder
-----------------------
* ``test_planted_truth.py`` asserts the facts the data generator planted
  (``planted_truth.json``) are found -- the outputs are RIGHT, where this file
  proves they did not CHANGE.  Both read the same run (``conftest.pipeline_run``).
* ``diff_runs.py`` compares any two output trees you made yourself (two
  branches, two commits) with the same comparator and contract.

``--update-golden`` snapshots the CURRENT branch's outputs as the new golden;
use it after an intentional output change and commit the golden with the code.
See ``tests/regression/README.md``.

Two scenarios are parametrized over via ``conftest.SCENARIOS``:

  * ``longitudinal`` — pre/post timepoint design, ``example_data_longitudinal/``
  * ``single``       — single-timepoint design, ``example_data_single/``

Both run in **profiles mode** (``config_generated_<dt>.yml``): the bundled
profiles carry the planted facts recorded in ``planted_truth.json``, which
the BAM-mode path (``config_with_bams_<dt>.yml``, the tutorial's Option B)
cannot, because its BAMs are simulated from the unmodified reference.

Run in place
------------
The test runs the pipeline **directly inside the committed
``example_data_<scenario>/`` directory** (no tmp staging).  This matches
the generator's design choice of writing absolute paths into the configs,
so the pipeline writes outputs at the location ``output.root_dir``
specifies — which is ``<example_data_dir>/output_bams/<data_type>/``.

The test wipes ``output_bams/`` and ``.snakemake/`` from the example dir
before each run so consecutive runs always start clean.  Both directories
are gitignored, so the dirty committed dir is transient.

Run with:

    # default — longitudinal only
    pytest tests/regression/ -v

    # single only
    pytest tests/regression/ --regression-scenario single -v

    # both
    pytest tests/regression/ --regression-scenario all -v

Regenerate the snapshot after an intentional change:

    pytest tests/regression/ --update-golden --regression-scenario all -v

The test is marked ``regression`` so it can be deselected from the fast
unit-test loop:

    pytest tests/ -m "not regression"
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from typing import Iterable

import pytest

from tests.regression.compare import compare_trees, forbidden_present
from tests.regression.conftest import SCENARIOS

logger = logging.getLogger(__name__)


# Files we intentionally do NOT track in the snapshot — Snakemake bookkeeping
# and per-run metadata that varies even with deterministic inputs.
_EXCLUDED_SUFFIXES = (".log", ".snakemake_timestamp")
_EXCLUDED_PREFIXES = ("config_runtime_",)
# .ffn = the dN/dS step's reconstructed per-subject gene sequences (text).
_TRACKED_SUFFIXES = (".tsv", ".tsv.gz", ".parquet", ".json", ".ffn")


def _copy_golden_payload(
    src: Path, dst: Path, exclude_substrings: Iterable[str] = ()
) -> None:
    """Copy the subset of pipeline outputs we track in the snapshot.

    ``exclude_substrings`` is the scenario's comparison carve-out (the cache
    parquet files, the profile copies); applying it here too keeps those
    files out of the committed golden instead of merely skipping them.
    """
    excl = tuple(exclude_substrings)
    for f in src.rglob("*"):
        if not f.is_file():
            continue
        if f.name.startswith(_EXCLUDED_PREFIXES):
            continue
        if f.name.endswith(_EXCLUDED_SUFFIXES):
            continue
        if not f.name.endswith(_TRACKED_SUFFIXES):
            continue
        rel = f.relative_to(src)
        if excl and any(s in rel.as_posix() for s in excl):
            continue
        target = dst / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(f, target)


@pytest.mark.regression
def test_pipeline_matches_golden(
    scenario,
    pipeline_run,
    golden_root,
    update_golden,
):
    """Run one scenario of the example pipeline and diff against its golden.

    Two modes, picked via ``--update-golden``:

    * **Normal mode** (default): runs the pipeline on the CURRENT branch
      using the committed
      ``docs/source/examples/example_data_<scenario>/`` directory, then
      diffs its outputs against the committed golden snapshot at
      ``tests/regression/golden/<scenario>/``.
    * **Update mode** (``--update-golden``): snapshots the CURRENT branch's
      outputs into the golden directory.  Use this when your branch
      intentionally changes pipeline outputs and you want those outputs to
      become the new baseline, and commit the golden with the change.

    The pipeline must be installed in the active conda environment
    (``pip install -e .`` from the repo root, which is what the alleleflux
    env already provides).  The run is the ``pipeline_run`` fixture
    (``conftest.py``): one pipeline run per scenario per session, shared
    with ``test_planted_truth.py``.
    """
    bundle = SCENARIOS[scenario]
    # The run itself lives in conftest.pipeline_run (one per scenario per
    # session) so test_planted_truth.py can read the same outputs.
    pipeline_data_dir = pipeline_run.pipeline_data_dir

    # Trees this scenario must NOT produce (e.g. the longitudinal-only
    # dN/dS / baseline-presence / strain-turnover outputs on single data).
    # Checked before the golden diff and in update mode too: a golden must
    # never be captured from a run that wrote a forbidden tree.
    forbidden = forbidden_present(pipeline_data_dir, bundle.get("expected_absent", ()))
    if forbidden:
        pytest.fail(
            f"[scenario={scenario}] outputs that must not exist were written:\n"
            + "\n".join(forbidden)
        )

    golden_dir = golden_root / bundle["golden_subdir"]

    # --update-golden mode: this is NOT a check, it is a snapshot refresh
    # FROM the current branch.  Blow away the old golden, repopulate from
    # this run, and skip (so the run is not also reported as a pass/fail
    # against itself).
    if update_golden:
        if golden_dir.exists():
            shutil.rmtree(golden_dir)
        golden_dir.mkdir(parents=True)
        _copy_golden_payload(
            pipeline_data_dir, golden_dir, bundle.get("exclude_substrings", ())
        )
        pytest.skip(
            f"Golden snapshot refreshed at {golden_dir} from the CURRENT branch; "
            f"review `git diff --stat tests/regression/golden` and commit it with the change."
        )

    assert golden_dir.exists(), (
        f"Golden snapshot missing for scenario '{scenario}' at {golden_dir}. "
        f"Run `pytest tests/regression/ --update-golden --regression-scenario {scenario}` to create it."
    )

    # The golden is captured from THIS branch, so the whole tracked tree is
    # compared at the comparator's strict defaults (FULL contract, _scope.py);
    # a scenario may still narrow the sweep or loosen the tolerance through
    # its bundle, which is how the retired main-derived mode used to work.
    tol = {k: bundle[k] for k in ("rtol", "atol") if k in bundle}
    failures = compare_trees(
        pipeline_data_dir,
        golden_dir,
        include_prefixes=bundle.get("compare_only"),
        allow_golden_only=bundle.get("allow_golden_only"),
        exclude_substrings=bundle.get("exclude_substrings"),
        **tol,
    )
    if failures:
        joined = "\n\n".join(failures)
        pytest.fail(
            f"[scenario={scenario}] {len(failures)} file(s) differ from the "
            f"golden snapshot:\n\n{joined}"
        )
