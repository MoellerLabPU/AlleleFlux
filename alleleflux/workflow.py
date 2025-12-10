#!/usr/bin/env python3
"""
AlleleFlux Workflow Orchestration Module.

This module handles the actual execution of Snakemake workflows,
including configuration validation, command building, and real-time
output streaming. It is separate from the CLI to maintain clean
separation of concerns.
"""

import logging
import multiprocessing
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import pkg_resources
import yaml

from alleleflux.scripts.utilities.logging_config import setup_logging

setup_logging()
logger = logging.getLogger(__name__)
# =============================================================================
# Path Helpers
# =============================================================================


def get_snakefile() -> Path:
    """
    Get the path to the packaged Snakefile.

    Returns:
        Path to the Snakefile bundled with the package.

    Raises:
        SystemExit: If the Snakefile cannot be found.
    """
    package_dir = Path(__file__).parent
    snakefile = package_dir / "smk_workflow" / "alleleflux_pipeline" / "Snakefile"

    if not snakefile.exists():
        logger.error(f"Snakefile not found at: {snakefile}")
        logger.error("This may indicate a broken installation.")
        sys.exit(1)

    return snakefile


def get_template_path() -> Path:
    """
    Get the path to the configuration template file.

    Returns:
        Path to the config.template.yml bundled with the package.
    """
    package_dir = Path(__file__).parent
    template = package_dir / "smk_workflow" / "config.template.yml"

    # Fallback to pkg_resources if Path-based lookup fails
    if not template.exists():
        try:
            return Path(
                pkg_resources.resource_filename(
                    "alleleflux", "smk_workflow/config.template.yml"
                )
            )
        except Exception:
            pass

    return template


# =============================================================================
# Configuration Validation
# =============================================================================


def validate_config(config_path: Path) -> dict:
    """
    Load and validate a configuration file.

    Args:
        config_path: Path to the configuration file.

    Returns:
        Parsed configuration dictionary.

    Raises:
        SystemExit: If configuration is invalid or missing required fields.
    """
    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        logger.error("Generate one with 'alleleflux init' or provide a valid path.")
        sys.exit(1)

    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        logger.error(f"Error parsing configuration file: {e}")
        sys.exit(1)

    # Validate required sections
    required_sections = [
        "input",
        "output",
        "analysis",
        "quality_control",
        "profiling",
        "statistics",
        "dnds",
    ]
    missing = [s for s in required_sections if s not in config]
    if missing:
        logger.error(f"Configuration missing required sections: {missing}")
        sys.exit(1)

    # Validate required input paths
    input_config = config.get("input", {})
    required_inputs = [
        "fasta_path",
        "prodigal_path",
        "metadata_path",
        "mag_mapping_path",
        "gtdb_path",
    ]
    missing_inputs = [i for i in required_inputs if not input_config.get(i)]
    if missing_inputs:
        logger.error(f"Configuration missing required input paths: {missing_inputs}")
        sys.exit(1)

    return config


# =============================================================================
# Snakemake Execution
# =============================================================================


def run_snakemake(cmd: str, logfile: Path = None, verbose: bool = False) -> int:
    """
    Run snakemake command with real-time output streaming.

    This follows the SnakePipes pattern of using subprocess.Popen with
    stdout streaming for immediate user feedback during long workflow runs.
    Handles keyboard interrupts gracefully by allowing snakemake to perform
    its own cleanup before exiting.

    Args:
        cmd: The snakemake command string to execute.
        logfile: Optional path to write command output to a log file.
        verbose: If True, print the command before execution.

    Returns:
        Exit code from snakemake (0 = success, 130 = interrupted).
    """
    if verbose:
        logger.info(f"Executing: {cmd}")

    # Open log file if specified
    log_handle = open(logfile, "w") if logfile else None

    try:
        # Run snakemake with real-time output streaming
        # Using shell=True following SnakePipes/Atlas pattern
        process = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )

        # Stream output line by line
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            if log_handle:
                log_handle.write(line)
                log_handle.flush()

        # Wait for process to complete
        process.wait()
        return process.returncode

    except KeyboardInterrupt:
        # Handle Ctrl+C gracefully - let snakemake clean up
        logger.warning("\nInterrupted by user. Waiting for snakemake to clean up...")
        if log_handle:
            log_handle.write("\nInterrupted by user (KeyboardInterrupt)\n")
            log_handle.flush()

        # Snakemake receives SIGINT automatically since it's in the same process group
        # Wait for it to finish its cleanup
        try:
            process.wait(timeout=30)  # Give snakemake 30 seconds to clean up
        except subprocess.TimeoutExpired:
            logger.warning("Snakemake cleanup timed out. Terminating...")
            process.terminate()
            process.wait()

        logger.info("Workflow interrupted. You can resume with the same command.")
        return 130  # Standard exit code for SIGINT

    finally:
        if log_handle:
            log_handle.close()


def build_snakemake_command(
    snakefile: Path,
    config_path: Path,
    working_dir: str,
    jobs: Optional[int] = None,
    profile: Optional[str] = None,
    dry_run: bool = False,
    extra_args: Optional[List[str]] = None,
) -> str:
    """
    Build the snakemake command string.

    Args:
        snakefile: Path to the Snakefile.
        config_path: Path to the configuration file.
        working_dir: Working directory for the workflow.
        jobs: Number of parallel jobs (None = auto-detect).
        profile: Path to Snakemake profile for cluster execution.
        dry_run: If True, add --dry-run flag.
        extra_args: Additional arguments to pass to snakemake.

    Returns:
        Complete snakemake command string.
    """
    cmd_parts = [
        "snakemake",
        f"--snakefile {snakefile}",
        f"--configfile {config_path}",
        f"--directory {working_dir}",
        "--rerun-incomplete",
        "--show-failed-logs",
    ]

    # Add jobs/cores
    if jobs is not None:
        cmd_parts.append(f"--jobs {jobs}")
    elif profile is None:
        # Local execution without profile: use all cores
        cmd_parts.append(f"--cores {multiprocessing.cpu_count()}")

    # Add profile for cluster execution
    if profile:
        cmd_parts.append(f"--profile {profile}")

    # Add dry run flag
    if dry_run:
        cmd_parts.append("--dry-run")

    # Add any extra snakemake arguments
    if extra_args:
        cmd_parts.extend(extra_args)

    return " ".join(cmd_parts)


def execute_workflow(
    config_file: str,
    working_dir: str = ".",
    jobs: Optional[int] = None,
    profile: Optional[str] = None,
    dry_run: bool = False,
    unlock: bool = False,
    verbose: bool = False,
    extra_args: Optional[List[str]] = None,
    version: str = "unknown",
) -> int:
    """
    Execute the AlleleFlux workflow.

    This is the main entry point for workflow execution, handling all
    setup, validation, and execution logic.

    Args:
        config_file: Path to the configuration file.
        working_dir: Working directory for the workflow.
        jobs: Number of parallel jobs.
        profile: Path to Snakemake profile for cluster execution.
        dry_run: If True, perform a dry run without executing jobs.
        unlock: If True, unlock the working directory and exit.
        verbose: If True, print the snakemake command before execution.
        extra_args: Additional arguments to pass to snakemake.
        version: AlleleFlux version string for logging.

    Returns:
        Exit code (0 = success).
    """
    # Get the Snakefile path
    snakefile = get_snakefile()

    # Validate configuration
    config_path = Path(config_file)
    config = validate_config(config_path)

    # Handle unlock
    if unlock:
        logger.info("Unlocking working directory...")
        unlock_cmd = f"snakemake --snakefile {snakefile} --unlock"
        return run_snakemake(unlock_cmd, verbose=verbose)

    # Get output directory from config for logging
    output_dir = Path(config.get("output", {}).get("root_dir", working_dir))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create log file with timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    logfile = output_dir / f"alleleflux_{timestamp}.log"

    # Build snakemake command
    cmd = build_snakemake_command(
        snakefile=snakefile,
        config_path=config_path,
        working_dir=working_dir,
        jobs=jobs,
        profile=profile,
        dry_run=dry_run,
        extra_args=extra_args,
    )

    # Log startup info
    logger.info("=" * 60)
    logger.info("AlleleFlux Workflow Starting")
    logger.info("=" * 60)
    logger.info(f"Version: {version}")
    logger.info(f"Snakefile: {snakefile}")
    logger.info(f"Config: {config_path}")
    logger.info(f"Working directory: {working_dir}")
    if profile:
        logger.info(f"Profile: {profile}")
    logger.info(f"Log file: {logfile}")
    if dry_run:
        logger.warning("DRY RUN MODE - No jobs will be executed")
    logger.info("=" * 60)

    # Write command to log
    with open(logfile, "w") as f:
        f.write(f"AlleleFlux run started at {timestamp}\n")
        f.write(f"Snakemake command: {cmd}\n\n")

    # Run snakemake
    try:
        exit_code = run_snakemake(cmd, logfile=logfile, verbose=verbose)
    except KeyboardInterrupt:
        # This catches any keyboard interrupt not caught in run_snakemake
        logger.warning("\nWorkflow interrupted by user.")
        return 130

    # Report result
    if exit_code == 0:
        logger.info("=" * 60)
        logger.info("AlleleFlux workflow completed successfully!")
        logger.info(f"Output directory: {output_dir}")
        logger.info(f"Log file: {logfile}")
        logger.info("=" * 60)
    elif exit_code == 130:
        # Already logged in run_snakemake, just return
        pass
    else:
        logger.error("=" * 60)
        logger.error(f"AlleleFlux workflow failed with exit code {exit_code}")
        logger.error(f"Check log file for details: {logfile}")
        logger.error("=" * 60)

    return exit_code
