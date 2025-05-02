#!/usr/bin/env python3
# AlleleFlux Pipeline Runner
import argparse
import logging
import multiprocessing as mp
import os
import signal
import subprocess
import sys
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, Dict, Optional

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


def run_snakemake(
    snakefile: str,
    config_file: str,
    threads: int,
    profile: Optional[str] = None,
    dry_run: bool = False,
) -> bool:
    """Run Snakemake with proper process group isolation."""
    global child_proc
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

    logger.info(f"Running command: {' '.join(cmd)}")

    try:
        # Start Snakemake in a new process group with stdout/stderr passed through
        child_proc = subprocess.Popen(
            cmd,
            preexec_fn=os.setsid,  # Isolate process group
            stdout=None,  # Use parent's stdout (pass through)
            stderr=None,  # Use parent's stderr (pass through)
            bufsize=1,  # Line buffered
            universal_newlines=True,  # Text mode
        )
        child_proc.wait()
        return child_proc.returncode == 0
    except Exception as e:
        logger.error(f"Exception running {snakefile}: {e}")
        return False
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
        "-v",
        "--version",
        action="version",
        version=f"AlleleFlux v{version("AlleleFlux")}",
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
    if args.min_sample_num is not None:
        config["quality_control"]["min_sample_num"] = args.min_sample_num
    if args.breadth_threshold is not None:
        config["quality_control"]["breadth_threshold"] = args.breadth_threshold

    return config


def write_runtime_config(config: Dict[str, Any], outDir) -> str:
    """Write updated config to a runtime file and return its path."""
    timestamp = datetime.now().strftime("%m-%d-%Y_%H-%M-%S")
    runtime_config = f"alleleflux_config_{timestamp}.yml"
    runtime_config_path = os.path.join(outDir, runtime_config)

    try:
        with open(runtime_config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        logger.info(f"Created runtime config file: {runtime_config_path}")
        return runtime_config_path
    except Exception as e:
        logger.error(f"Error writing runtime configuration file: {e}")
        sys.exit(1)


def main():
    # Register signal handlers for graceful termination
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    # Parse command-line arguments
    args = parse_arguments()

    # Load configuration from specified file
    logger.info(f"Reading config file from: {args.config}")
    config = load_config(args.config)

    # Update configuration with command-line arguments
    updated_config = update_config(config, args)

    # Check for required parameters
    missing_params = []
    for param_name, param_value in [
        ("BAM directory", updated_config.get("input", {}).get("bam_dir")),
        ("FASTA file", updated_config.get("input", {}).get("fasta_path")),
        ("Prodigal file", updated_config.get("input", {}).get("prodigal_path")),
        ("Metadata file", updated_config.get("input", {}).get("metadata_path")),
        ("Output directory", updated_config.get("output", {}).get("root_dir")),
    ]:
        if not param_value:
            missing_params.append(param_name)

    if missing_params:
        logger.error("Missing required parameters: " + ", ".join(missing_params))
        logger.error(
            "Please specify these via command-line arguments or in the config file."
        )
        sys.exit(1)

    # Write runtime config file
    runtime_config_path = write_runtime_config(
        updated_config, updated_config.get("output", {}).get("root_dir")
    )
    logger.info(f"Created runtime configuration file: {runtime_config_path}")

    # Workflow 1: Step 1 (Preprocessing)
    logger.info("=" * 80)
    logger.info("Starting Step 1: Sample Profiling and Eligibility Table Generation")
    logger.info("=" * 80)

    if not run_snakemake(
        "smk_workflow/step1.smk",
        runtime_config_path,
        args.threads,
        args.profile,
        args.dry_run,
    ):
        logger.error("Step 1 failed! Aborting.")
        sys.exit(1)

    # Workflow 2: Step 2 (Allele Analysis)
    logger.info("=" * 80)
    logger.info("Starting Step 2: Allele Analysis and Scoring")
    logger.info("=" * 80)

    if not run_snakemake(
        "smk_workflow/step2.smk",
        runtime_config_path,
        args.threads,
        args.profile,
        args.dry_run,
    ):
        logger.error("Step 2 failed!")
        sys.exit(1)

    logger.info("Pipeline completed successfully!")
    logger.info(f"Runtime configuration preserved at: {runtime_config_path}")


if __name__ == "__main__":
    main()
