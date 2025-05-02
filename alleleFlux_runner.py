#!/usr/bin/env python3
# AlleleFlux Pipeline Runner
# A robust script to run the two-step Snakemake workflow with proper interruption handling

import argparse
import logging
import multiprocessing as mp
import os
import re
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import yaml

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("AlleleFlux")

# Global variables
child_proc = None  # Global variable to track the Snakemake process
CONFIG_PATH = os.path.join("smk_workflow", "config.yml")
STEP1_SNAKEFILE = os.path.join("smk_workflow", "step1.smk")
STEP2_SNAKEFILE = os.path.join("smk_workflow", "step2.smk")
INCOMPLETE_FILES_PATTERN = re.compile(
    r"Incomplete files:\s*(.*?)(?:\n\n|\Z)", re.DOTALL
)


def signal_handler(sig, frame):
    """Terminate Snakemake and its entire process group."""
    global child_proc
    logger.info(
        f"Received {signal.Signals(sig).name}, terminating Snakemake and subprocesses..."
    )

    if child_proc is not None:
        try:
            # Terminate the entire process group
            pgid = os.getpgid(child_proc.pid)
            os.killpg(pgid, sig)
            # Wait for process group to terminate
            child_proc.wait(timeout=15)
            logger.info("Snakemake and subprocesses terminated.")
        except subprocess.TimeoutExpired:
            logger.warning("Forcing immediate termination with SIGKILL...")
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            logger.info("Process group already terminated.")
        except Exception as e:
            logger.error(f"Error terminating process group: {e}")
    sys.exit(1)


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    try:
        with open(config_path, "r") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)
    except yaml.YAMLError as e:
        logger.error(f"Error parsing configuration file: {e}")
        sys.exit(1)


def validate_path(path: str, is_dir: bool = False, should_exist: bool = True) -> bool:
    """Validate if a path exists and is of the correct type."""
    path_obj = Path(path)

    if should_exist:
        if not path_obj.exists():
            logger.warning(f"Path does not exist: {path}")
            return False

        if is_dir and not path_obj.is_dir():
            logger.warning(f"Path is not a directory: {path}")
            return False
        elif not is_dir and not path_obj.is_file():
            logger.warning(f"Path is not a file: {path}")
            return False

    return True


def validate_config(config: Dict[str, Any]) -> bool:
    """Validate configuration parameters."""
    # Check required input paths
    required_inputs = [
        ("bam_dir", True),
        ("fasta_path", False),
        ("prodigal_path", False),
        ("metadata_path", False),
    ]

    validation_passed = True

    for input_name, is_dir in required_inputs:
        path = config.get("input", {}).get(input_name)
        if not path:
            logger.error(f"Missing required input: {input_name}")
            validation_passed = False
            continue

        # Just warn but don't fail validation if paths don't exist
        validate_path(path, is_dir, should_exist=True)

    # Check output directory
    output_dir = config.get("output", {}).get("root_dir")
    if not output_dir:
        logger.error("Missing output directory configuration")
        validation_passed = False

    # Validate data_type
    data_type = config.get("analysis", {}).get("data_type")
    if data_type not in ["single", "longitudinal"]:
        logger.error(
            f"Invalid data_type: {data_type}. Must be 'single' or 'longitudinal'"
        )
        validation_passed = False

    return validation_passed


def extract_incomplete_files(error_output: str) -> List[str]:
    """Extract incomplete file paths from Snakemake error output."""
    match = INCOMPLETE_FILES_PATTERN.search(error_output)
    if match:
        files_str = match.group(1).strip()
        return [f.strip() for f in files_str.split("\n")]
    return []


def run_snakemake(
    snakefile: str,
    config_file: str,
    threads: int,
    profile: Optional[str] = None,
    dry_run: bool = False,
    rerun_incomplete: bool = True,  # Default to True
    force_all: bool = False,
) -> Tuple[bool, List[str]]:
    """Run Snakemake with proper process group isolation and error handling."""
    global child_proc

    # Build the command
    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--configfile",
        config_file,
        "--cores",
        str(threads),
    ]

    # Add optional parameters
    if profile:
        cmd.extend(["--profile", profile])

    if dry_run:
        cmd.append("-n")
        logger.info(f"{' DRY RUN ':-^80}")
    else:
        logger.info(f"{' RUNNING ':-^80}")

    if rerun_incomplete:
        cmd.append("--rerun-incomplete")
        logger.info("Running with --rerun-incomplete flag")

    if force_all:
        cmd.append("--forceall")
        logger.info("Running with --forceall flag")

    logger.info(f"Running command: {' '.join(cmd)}")

    incomplete_files = []

    try:
        # Start Snakemake in a new process group with direct stdout/stderr
        process = subprocess.Popen(
            cmd,
            preexec_fn=os.setsid,  # Isolate process group
            stdout=None,  # Pass through to parent stdout
            stderr=None,  # Pass through to parent stderr
            bufsize=1,  # Line buffered
            universal_newlines=True,  # Text mode
        )

        child_proc = process
        process.wait()  # Wait for completion

        if process.returncode == 0:
            logger.info("Snakemake completed successfully")
            return True, []

        logger.error(f"Snakemake failed with return code {process.returncode}")
        return False, incomplete_files

    except Exception as e:
        logger.error(f"Exception running {snakefile}: {e}")
        return False, incomplete_files
    finally:
        child_proc = None  # Reset after completion


def parse_arguments():
    """Parse command-line arguments."""
    # Load config file if it exists
    try:
        config = load_config(CONFIG_PATH)
    except (FileNotFoundError, yaml.YAMLError):
        # If config file doesn't exist or has errors, use empty dictionary
        config = {"input": {}, "output": {}, "analysis": {}, "quality_control": {}}

    parser = argparse.ArgumentParser(
        description="AlleleFlux: A tool for fine grained evolutionary analysis of microbial populations and communities.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Basic arguments
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Perform a dry run (don't execute commands)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"AlleleFlux Runner v1.0.3",
        help="Show version information and exit",
    )
    runtime_group = parser.add_argument_group("Runtime configuration")

    # Runtime configuration
    runtime_group.add_argument(
        "--config", default=CONFIG_PATH, help="Path to configuration file"
    )
    runtime_group.add_argument(
        "--profile", default=None, help="Path to Snakemake profile directory"
    )
    runtime_group.add_argument(
        "--threads",
        type=int,
        default=mp.cpu_count(),
        help="Number of CPU threads to use",
    )
    runtime_group.add_argument(
        "--force-all",
        action="store_true",
        help="Force rerunning all steps (equivalent to snakemake --forceall)",
    )
    runtime_group.add_argument(
        "--no-rerun-incomplete",
        action="store_true",
        help="Disable automatic rerunning of jobs with incomplete output files",
    )
    runtime_group.add_argument(
        "--step1-only",
        action="store_true",
        help="Run only step 1 of the workflow",
    )
    runtime_group.add_argument(
        "--step2-only",
        action="store_true",
        help="Run only step 2 of the workflow",
    )

    # Input files
    input_group = parser.add_argument_group("Input Files")
    input_group.add_argument(
        "--bam-dir",
        default=config.get("input", {}).get("bam_dir"),
        help="Directory containing BAM files",
    )
    input_group.add_argument(
        "--fasta-path",
        default=config.get("input", {}).get("fasta_path"),
        help="Reference fasta file",
    )
    input_group.add_argument(
        "--prodigal-path",
        default=config.get("input", {}).get("prodigal_path"),
        help="Prodigal output file",
    )
    input_group.add_argument(
        "--metadata-path",
        default=config.get("input", {}).get("metadata_path"),
        help="Sample metadata file",
    )

    # Output directory
    output_group = parser.add_argument_group("Output Files")
    output_group.add_argument(
        "--output-dir",
        default=config.get("output", {}).get("root_dir"),
        help="Root directory for output files",
    )

    # Analysis parameters
    analysis_group = parser.add_argument_group("Analysis Parameters")
    analysis_group.add_argument(
        "--data-type",
        default=config.get("analysis", {}).get("data_type"),
        choices=["single", "longitudinal"],
        help="Data type for analysis approach",
    )
    analysis_group.add_argument(
        "--use-lmm",
        action="store_true",
        default=config.get("analysis", {}).get("use_lmm", True),
        help="Enable Linear Mixed Models analysis",
    )
    analysis_group.add_argument(
        "--use-significance-tests",
        action="store_true",
        default=config.get("analysis", {}).get("use_significance_tests", True),
        help="Enable significance tests",
    )

    # Quality control parameters
    qc_group = parser.add_argument_group("Quality Control Parameters")
    qc_group.add_argument(
        "--min-sample-num",
        type=int,
        default=config.get("quality_control", {}).get("min_sample_num"),
        help="Minimum number of replicates required to calculate p-values",
    )
    qc_group.add_argument(
        "--breadth-threshold",
        type=float,
        default=config.get("quality_control", {}).get("breadth_threshold"),
        help="Minimum breadth to analyze a MAG",
    )

    return parser.parse_args()


def update_config(config: Dict[str, Any], args) -> Dict[str, Any]:
    """Update configuration with command-line arguments."""
    # Ensure input section exists
    if "input" not in config:
        config["input"] = {}

    # Update input paths if provided on command line
    if args.bam_dir:
        config["input"]["bam_dir"] = args.bam_dir
    if args.fasta_path:
        config["input"]["fasta_path"] = args.fasta_path
    if args.prodigal_path:
        config["input"]["prodigal_path"] = args.prodigal_path
    if args.metadata_path:
        config["input"]["metadata_path"] = args.metadata_path

    # Ensure output section exists
    if "output" not in config:
        config["output"] = {}

    # Update output directory if provided
    if args.output_dir:
        config["output"]["root_dir"] = args.output_dir

    # Ensure analysis section exists
    if "analysis" not in config:
        config["analysis"] = {}

    # Update analysis parameters if provided
    if args.data_type:
        config["analysis"]["data_type"] = args.data_type
    config["analysis"]["use_lmm"] = args.use_lmm
    config["analysis"]["use_significance_tests"] = args.use_significance_tests

    # Ensure quality_control section exists
    if "quality_control" not in config:
        config["quality_control"] = {}

    # Update quality control parameters if provided
    if args.min_sample_num:
        config["quality_control"]["min_sample_num"] = args.min_sample_num
    if args.breadth_threshold:
        config["quality_control"]["breadth_threshold"] = args.breadth_threshold

    return config


def save_config(config: Dict[str, Any], config_path: str) -> bool:
    """Save configuration to YAML file."""
    try:
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        return True
    except Exception as e:
        logger.error(f"Error saving configuration: {e}")
        return False


def write_runtime_config(config: Dict[str, Any], out_dir: str) -> str:
    """Write updated config to a runtime file and return its path."""
    timestamp = datetime.now().strftime("%m-%d-%Y_%H:%M:%S")
    runtime_config = f"alleleflux_config_{timestamp}.yml"
    runtime_config_path = os.path.join(out_dir, runtime_config)

    try:
        os.makedirs(out_dir, exist_ok=True)
        with open(runtime_config_path, "w") as f:
            # Use flow style for lists to maintain original format
            yaml.dump(
                config, f, default_flow_style=None, sort_keys=False, default_style=None
            )
        logger.info(f"Created runtime config file: {runtime_config_path}")
        return runtime_config_path
    except Exception as e:
        logger.error(f"Error writing runtime configuration file: {e}")
        sys.exit(1)


def main():
    """Main entry point."""
    # Register signal handlers
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    # Parse command-line arguments
    args = parse_arguments()

    # Load configuration
    config = load_config(args.config)

    # Update configuration with command-line arguments
    config = update_config(config, args)

    # Validate configuration
    if not validate_config(config):
        logger.error("Configuration validation failed. Please check your settings.")
        sys.exit(1)

    # Create runtime config file in output directory
    output_dir = config.get("output", {}).get("root_dir")
    if not output_dir:
        logger.error("Output directory not specified in configuration")
        sys.exit(1)

    runtime_config_path = write_runtime_config(config, output_dir)

    # Run the workflow
    success = True

    # Determine which steps to run
    run_step1 = not args.step2_only
    run_step2 = not args.step1_only

    # Determine whether to use --rerun-incomplete flag
    rerun_incomplete = not args.no_rerun_incomplete

    # Run Step 1: Sample Profiling and Eligibility Table Generation
    if run_step1:
        logger.info(
            "Starting Step 1: Sample Profiling and Eligibility Table Generation"
        )
        step1_success, incomplete_files = run_snakemake(
            STEP1_SNAKEFILE,
            runtime_config_path,
            args.threads,
            args.profile,
            args.dry_run,
            rerun_incomplete,
            args.force_all,
        )

        if not step1_success:
            if incomplete_files:
                logger.error(
                    f"Step 1 failed with incomplete files: {', '.join(incomplete_files)}"
                )
            else:
                logger.error("Step 1 failed")

            # Only exit if step 2 is not requested to run independently
            if not args.step2_only:
                sys.exit(1)
        else:
            logger.info("Step 1 completed successfully")

        success = success and step1_success

    # Run Step 2: Allele Analysis and Scoring
    if run_step2:
        logger.info("Starting Step 2: Allele Analysis and Scoring")
        step2_success, incomplete_files = run_snakemake(
            STEP2_SNAKEFILE,
            runtime_config_path,
            args.threads,
            args.profile,
            args.dry_run,
            rerun_incomplete,
            args.force_all,
        )

        if not step2_success:
            if incomplete_files:
                logger.error(
                    f"Step 2 failed with incomplete files: {', '.join(incomplete_files)}"
                )
            else:
                logger.error("Step 2 failed")

            success = False
        else:
            logger.info("Step 2 completed successfully")

        success = success and step2_success

    # Final status
    if success:
        logger.info("AlleleFlux workflow completed successfully")
    else:
        logger.error("AlleleFlux workflow completed with errors")
        sys.exit(1)


if __name__ == "__main__":
    main()
