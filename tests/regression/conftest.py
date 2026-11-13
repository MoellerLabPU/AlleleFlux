"""
Pytest configuration for the end-to-end regression test.

Adds three custom CLI options:

  --update-golden        Re-run the pipeline on the example dataset and overwrite
                         the golden snapshot in tests/regression/golden/. Use after
                         an intentional change to pipeline outputs.

  --regression-scenario  Which scenario(s) to run.  One of
                         ``longitudinal`` (default), ``single``, or ``all``.
                         Each scenario is a (data type, dataset directory,
                         config, golden subdir) bundle — see ``SCENARIOS``
                         below.

  --regression-config    Optional explicit override for the config path used in
                         the run.  Rare — primarily for ad-hoc comparisons.
                         When given, it replaces the scenario's default config.

The ``SCENARIOS`` mapping is the single source of truth for what the test
exercises.  Add a new entry to introduce a new scenario; ``test_pipeline_e2e``
parametrizes over the entries selected by ``--regression-scenario``.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest
import yaml

# The comparison contract, shared with diff_runs.py so both apply the identical
# rules.  See _scope.py.
from tests.regression._scope import (
    EXPECTED_ABSENT as _EXPECTED_ABSENT,
    FULL_EXCLUDE_SUBSTRINGS as _FULL_EXCLUDE_SUBSTRINGS,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
GOLDEN_ROOT = REPO_ROOT / "tests/regression/golden"
EXAMPLES_DIR = REPO_ROOT / "docs/source/examples"


# Single source of truth for which scenarios the regression test covers.
# Each entry resolves to a fully-qualified bundle: where the inputs live,
# which config drives the run, where the pipeline writes its outputs, and
# where the committed golden snapshot lives.
#
# Both scenarios run in PROFILES mode (config_generated_<dt>.yml): the
# bundled profiles carry the planted facts (planted_truth.json), while the
# BAMs are simulated from the unmodified reference and carry none of them.
# The golden is captured from this branch, so the FULL contract applies:
# every tracked file in the output tree, strict tolerance, minus the
# exclusions in _scope.FULL_EXCLUDE_SUBSTRINGS; plus the per-scenario list of
# trees that must NOT exist (_scope.EXPECTED_ABSENT).
SCENARIOS: dict[str, dict] = {
    "longitudinal": {
        "data_type": "longitudinal",
        "example_data_dir": str(EXAMPLES_DIR / "example_data_longitudinal"),
        "config_name": "config_generated_longitudinal.yml",
        "pipeline_data_dir": "longitudinal",
        "golden_subdir": "longitudinal",
        "exclude_substrings": _FULL_EXCLUDE_SUBSTRINGS,
        "expected_absent": _EXPECTED_ABSENT["longitudinal"],
    },
    "single": {
        "data_type": "single",
        "example_data_dir": str(EXAMPLES_DIR / "example_data_single"),
        "config_name": "config_generated_single.yml",
        # The pipeline writes the single data_type under "single_timepoint"
        # (common.smk: DATA_TYPE=="single" -> OUTDIR/single_timepoint), NOT
        # "single".  Must match or the test can't find the outputs.
        "pipeline_data_dir": "single_timepoint",
        "golden_subdir": "single",
        "exclude_substrings": _FULL_EXCLUDE_SUBSTRINGS,
        "expected_absent": _EXPECTED_ABSENT["single"],
    },
}


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Overwrite the regression golden snapshot with the pipeline's "
        "current outputs.  Only use after an intentional output change.",
    )
    parser.addoption(
        "--regression-scenario",
        action="store",
        default="longitudinal",
        choices=[*SCENARIOS.keys(), "all"],
        help="Which scenario(s) to run.  ``all`` runs every entry in SCENARIOS.",
    )
    parser.addoption(
        "--regression-config",
        action="store",
        default=None,
        help="Optional override of the config path for the chosen scenario.  "
        "Defaults to the config bundled in the scenario's example_data_* dir.",
    )


def _selected_scenarios(config) -> list[str]:
    """Resolve --regression-scenario to a list of scenario names."""
    val = config.getoption("--regression-scenario")
    if val == "all":
        return list(SCENARIOS.keys())
    return [val]


def pytest_generate_tests(metafunc):
    """Parametrize the regression test over the requested scenarios.

    Tests that take a ``scenario`` parameter get one invocation per selected
    scenario.  Other tests are unaffected.
    """
    if "scenario" in metafunc.fixturenames:
        scenarios = _selected_scenarios(metafunc.config)
        metafunc.parametrize("scenario", scenarios, ids=scenarios)


@pytest.fixture(scope="session")
def update_golden(request) -> bool:
    return bool(request.config.getoption("--update-golden"))


@pytest.fixture
def regression_config_override(request) -> Path | None:
    """Optional explicit config-path override.  None when not set."""
    val = request.config.getoption("--regression-config")
    return Path(val).resolve() if val else None


@pytest.fixture(scope="session")
def golden_root() -> Path:
    return GOLDEN_ROOT


# =============================================================================
# One pipeline run per scenario, shared by every test in this folder
# =============================================================================

logger = logging.getLogger(__name__)

# Directory names the pipeline writes inside the example_data_* dir that we
# must wipe before each run (Snakemake reuses outputs based on timestamps,
# so a leftover from a previous run could silently masquerade as a fresh one).
_TRANSIENT_OUTPUT_DIRS = (
    "output_bams",  # written by config_with_bams_<dt>.yml
    "output",  # written by config_generated_<dt>.yml (profiles mode)
    ".snakemake",  # Snakemake's own bookkeeping
)


@dataclass(frozen=True)
class PipelineRun:
    """Where one scenario's run read its inputs and wrote its outputs."""

    scenario: str
    example_dir: Path  # the dataset folder; cwd of the run, holds planted_truth.json
    config_path: Path
    output_root: Path  # config output.root_dir, resolved
    pipeline_data_dir: Path  # output_root / <longitudinal|single_timepoint>


def _resolve_config_path(scenario: str, override: Path | None) -> Path:
    """The ``--regression-config`` override when given, else the scenario's bundled config."""
    if override is not None:
        return override
    bundle = SCENARIOS[scenario]
    return Path(bundle["example_data_dir"]) / bundle["config_name"]


def _resolve_output_root(config_path: Path) -> Path:
    """``output.root_dir`` from the config, resolved against the config's folder when relative."""
    with config_path.open() as f:
        cfg = yaml.safe_load(f)
    root_dir = Path(cfg["output"]["root_dir"])
    if not root_dir.is_absolute():
        root_dir = (config_path.parent / root_dir).resolve()
    return root_dir


def _wipe_transient_outputs(example_dir: Path, output_root: Path) -> None:
    """Remove anything the pipeline might have left from a previous run."""
    if output_root.exists():
        shutil.rmtree(output_root)
    for name in _TRANSIENT_OUTPUT_DIRS:
        p = example_dir / name
        if p.exists():
            shutil.rmtree(p)


def run_scenario(scenario: str, config_override: Path | None = None) -> PipelineRun:
    """Run one scenario's pipeline from a clean slate and return where it wrote.

    Runs ``alleleflux run`` with cwd = the dataset folder (the configs and
    metadata use paths relative to it) and 4 threads, after wiping every
    output a previous run may have left.  Fails the calling test with the
    pipeline's stdout/stderr tail when the run exits non-zero.
    """
    bundle = SCENARIOS[scenario]
    config_path = _resolve_config_path(scenario, config_override)
    assert config_path.exists(), f"Regression config not found: {config_path}"
    example_dir = config_path.parent
    output_root = _resolve_output_root(config_path)
    pipeline_data_dir = output_root / bundle["pipeline_data_dir"]
    _wipe_transient_outputs(example_dir, output_root)

    cmd = ["alleleflux", "run", "--config", str(config_path), "--threads", "4"]
    logger.info(
        f"Running pipeline for scenario={scenario}: {' '.join(cmd)} (cwd={example_dir})"
    )
    result = subprocess.run(cmd, cwd=example_dir, capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(
            f"Pipeline failed (scenario={scenario}):\n"
            f"  exit code: {result.returncode}\n"
            f"  stdout (tail):\n{result.stdout[-2000:]}\n"
            f"  stderr (tail):\n{result.stderr[-2000:]}"
        )
    assert pipeline_data_dir.exists(), (
        f"Pipeline data dir missing for scenario '{scenario}': {pipeline_data_dir}\n"
        f"Expected the pipeline to write outputs under {output_root}/."
    )
    return PipelineRun(
        scenario, example_dir, config_path, output_root, pipeline_data_dir
    )


# (scenario, config override) -> the finished run.  A plain module dict rather
# than a session-scoped fixture because the fixture is parametrized through
# ``scenario`` (pytest_generate_tests), and a session fixture cannot take a
# function-scoped parameter; the dict gives the same "run once per scenario"
# with no scope juggling.
_RUNS: dict[tuple[str, Path | None], PipelineRun] = {}


@pytest.fixture
def pipeline_run(scenario, regression_config_override) -> PipelineRun:
    """The scenario's pipeline run, executed once per session and then reused.

    Both ``test_pipeline_e2e.py`` (golden diff) and ``test_planted_truth.py``
    (planted-fact assertions) read the same outputs, so the ~4 minutes per
    scenario are paid once however many tests ask.
    """
    key = (scenario, regression_config_override)
    if key not in _RUNS:
        _RUNS[key] = run_scenario(scenario, regression_config_override)
    return _RUNS[key]
