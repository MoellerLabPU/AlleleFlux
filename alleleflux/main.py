#!/usr/bin/env python3
import argparse
import logging
import os
import signal
import subprocess
import sys
from datetime import datetime
from importlib import metadata, resources
from multiprocessing import cpu_count

import yaml

# --- Setup basic logging ---
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


# --- Custom YAML Formatting ---
# Create a custom list class to trigger flow-style YAML output for specific lists.
class FlowList(list):
    pass


def flow_list_representer(dumper, data):
    """Tells PyYAML to represent FlowList instances in flow style (e.g., [a, b])."""
    # By using 'tag:yaml.org,2002:seq', we are telling PyYAML to treat our FlowList object as a standard list for the purpose of formatting.
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)


# Register it globally
yaml.add_representer(FlowList, flow_list_representer)


def create_config(args):
    """Builds the configuration dictionary with command-line args overriding a base file."""
    config_data = {}
    if args.config:
        logging.info(f"Loading base configuration from: {args.config}")
        with open(args.config, "r") as f:
            config_data = yaml.safe_load(f)

    # Helper to safely set nested dictionary keys
    # Clean CLI override logic for nested YAML configs
    def set_nested(keys, value):
        d = config_data
        for key in keys[:-1]:
            d = d.setdefault(key, {})
        d[keys[-1]] = value

    # Set non-boolean command-line overrides
    if args.bam_dir:
        set_nested(["input", "bam_dir"], args.bam_dir)
    if args.fasta_path:
        set_nested(["input", "fasta_path"], args.fasta_path)
    if args.prodigal_path:
        set_nested(["input", "prodigal_path"], args.prodigal_path)
    if args.metadata_path:
        set_nested(["input", "metadata_path"], args.metadata_path)
    if args.gtdb_path:
        set_nested(["input", "gtdb_path"], args.gtdb_path)
    if args.mag_mapping_path:
        set_nested(["input", "mag_mapping_path"], args.mag_mapping_path)
    if args.root_dir:
        set_nested(["output", "root_dir"], args.root_dir)
    if args.data_type:
        set_nested(["analysis", "data_type"], args.data_type)
    if args.min_sample_num:
        set_nested(["quality_control", "min_sample_num"], args.min_sample_num)
    if args.breadth_threshold:
        set_nested(["quality_control", "breadth_threshold"], args.breadth_threshold)
    if args.p_value_threshold:
        set_nested(["statistics", "p_value_threshold"], args.p_value_threshold)
    if args.filter_type:
        set_nested(["statistics", "filter_type"], args.filter_type)
    if args.max_zero_count is not None:
        set_nested(["statistics", "max_zero_count"], args.max_zero_count)

    # --- Automatic Resource Allocation ---
    # Only perform this logic if the specific resource section is not already defined in a base config.
    if not config_data.get("resources", {}).get("cpus"):
        threads_for_parallel_jobs = min(args.cores, args.threads_per_job)
        logging.info(
            f"Setting default threads per parallelizable job to {threads_for_parallel_jobs}"
        )
        set_nested(["resources", "cpus", "threads_per_job"], threads_for_parallel_jobs)

    if not config_data.get("resources", {}).get("memory"):
        logging.info(
            f"Setting default memory tiers: low_mem={args.low_mem}MB, high_mem={args.high_mem}MB"
        )
        set_nested(["resources", "memory", "low_mem"], args.low_mem)
        set_nested(["resources", "memory", "high_mem"], args.high_mem)

    if not config_data.get("resources", {}).get("time"):
        logging.info(
            f"Setting default time tiers: low_time='{args.low_time}', high_time='{args.high_time}'"
        )
        set_nested(["resources", "time", "low_time"], args.low_time)
        set_nested(["resources", "time", "high_time"], args.high_time)

    # Layer defaults for booleans, then override with flags
    # The setdefault method only sets a value if the key does not already exist. So it will not overwrite if the key is already set in the config.
    config_data.setdefault("analysis", {}).setdefault("use_lmm", True)
    config_data.setdefault("analysis", {}).setdefault("use_significance_tests", True)
    config_data.setdefault("analysis", {}).setdefault("use_cmh", True)
    config_data.setdefault("analysis", {}).setdefault("allele_analysis_only", False)
    config_data.setdefault("statistics", {}).setdefault("preprocess_two_sample", True)
    config_data.setdefault("quality_control", {}).setdefault(
        "disable_zero_diff_filtering", False
    )

    if args.disable_lmm:
        config_data["analysis"]["use_lmm"] = False
    if args.disable_significance_tests:
        config_data["analysis"]["use_significance_tests"] = False
    if args.disable_cmh:
        config_data["analysis"]["use_cmh"] = False
    if args.allele_analysis_only:
        config_data["analysis"]["allele_analysis_only"] = True
    if args.disable_zero_diff_filtering:
        config_data["quality_control"]["disable_zero_diff_filtering"] = True
    if args.disable_preprocess_two_sample:
        config_data["statistics"]["preprocess_two_sample"] = False

    # Simplified logic for a single timepoint/group combination
    if args.timepoints:
        data_type = config_data.get("analysis", {}).get("data_type")
        if data_type == "longitudinal":
            combination = {"timepoint": FlowList(args.timepoints), "focus": args.focus}
        elif data_type == "single":
            combination = {"timepoint": FlowList(args.timepoints)}
        set_nested(["analysis", "timepoints_combinations"], [combination])

    if args.groups:
        set_nested(["analysis", "groups_combinations"], [FlowList(args.groups)])

    return config_data


def main():
    """Python entry point for the AlleleFlux pipeline."""
    parser = argparse.ArgumentParser(
        description="AlleleFlux: A tool to identify genomic targets of natural selection in bacterial communities.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # --- Grouped Arguments for Clarity ---
    control_group = parser.add_argument_group("Workflow Control")
    input_group = parser.add_argument_group("Input Files")
    output_group = parser.add_argument_group("Output Directory")
    analysis_group = parser.add_argument_group("Analysis Settings")
    qc_group = parser.add_argument_group("Quality Control & Statistics")
    resource_group = parser.add_argument_group("Resource Allocation")

    # Workflow Control
    version_str = f"AlleleFlux {metadata.version('AlleleFlux')}"
    control_group.add_argument("-v", "--version", action="version", version=version_str)
    control_group.add_argument(
        "-j",
        "--cores",
        type=int,
        default=cpu_count(),
        help="Number of CPU cores for Snakemake.",
    )
    control_group.add_argument(
        "--generate-config",
        type=str,
        metavar="FILE_PATH",
        help="Generate a config file based on other arguments and save it to FILE_PATH, then exit.",
    )
    control_group.add_argument(
        "-p",
        "--profile",
        type=str,
        help="Snakemake profile directory for cluster execution.",
    )
    control_group.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing jobs.",
    )
    control_group.add_argument(
        "--step1-only",
        action="store_true",
        help="Only run sample profiling and eligibility file generation (step 1).",
    )
    control_group.add_argument(
        "--step2-only",
        action="store_true",
        help="Only run allele analysis, significance testing and score calculation (step 2).",
    )
    control_group.add_argument(
        "-c", "--config", type=str, help="Path to a base YAML config file to override."
    )

    # Input Files
    input_group.add_argument(
        "--bam-dir", type=str, help="Directory containing sorted BAM files."
    )
    input_group.add_argument("--fasta-path", type=str, help="Reference FASTA file.")
    input_group.add_argument(
        "--prodigal-path",
        type=str,
        help="Prodigal nucleic acid output file for the reference FASTA.",
    )
    input_group.add_argument("--metadata-path", type=str, help="Path to metadata file.")
    input_group.add_argument("--gtdb-path", type=str, help="GTDB taxonomy file.")
    input_group.add_argument(
        "--mag-mapping-path", type=str, help="Path to MAG-to-contig mapping file."
    )

    # Output Directory
    output_group.add_argument(
        "--root-dir", type=str, help="Root directory for all output files."
    )

    # Analysis Settings
    analysis_group.add_argument(
        "--data-type",
        type=str,
        choices=["single", "longitudinal"],
        default="longitudinal",
        help="Type of analysis. This choice affects other arguments.",
    )
    analysis_group.add_argument(
        "--allele-analysis-only",
        action="store_true",
        default=False,
        help="If set, only runs allele analysis.",
    )
    analysis_group.add_argument(
        "--disable-lmm",
        action="store_true",
        default=False,
        help="Disable Linear Mixed Models analysis (default: enabled).",
    )
    analysis_group.add_argument(
        "--disable-significance-tests",
        action="store_true",
        default=False,
        help="Disable two-sample and single-sample tests (default: enabled).",
    )
    analysis_group.add_argument(
        "--disable-cmh",
        action="store_true",
        default=False,
        help="Disable CMH test for categorical data (default: enabled).",
    )

    analysis_group.add_argument(
        "--timepoints",
        nargs="+",
        help="Define the timepoint combination.\n(e.g., for single: --timepoints end)\n(e.g., for longitudinal: --timepoints pre post)",
    )
    analysis_group.add_argument(
        "--focus",
        type=str,
        help="Define the focus timepoint for a longitudinal comparison.\n(Required for longitudinal mode)",
    )
    analysis_group.add_argument(
        "--groups",
        nargs="+",
        help="Define the group comparison.\n(e.g., --groups fat control)",
    )

    # Quality Control & Statistics
    qc_group.add_argument(
        "--min-sample-num",
        type=int,
        default=4,
        help="Minimum number of samples required.",
    )
    qc_group.add_argument(
        "--breadth-threshold",
        type=float,
        default=0.1,
        help="Minimum coverage breadth across genome.",
    )
    qc_group.add_argument(
        "--p-value-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for significance.",
    )
    qc_group.add_argument(
        "--filter-type",
        type=str,
        default="t-test",
        choices=["t-test", "wilcoxon", "both", "either"],
        help="Filter type for two-sample preprocessing.",
    )
    qc_group.add_argument(
        "--max-zero-count",
        type=int,
        default=4,
        help="Maximum number of zero values allowed in single-sample test.",
    )
    qc_group.add_argument(
        "--disable-zero-diff-filtering",
        action="store_true",
        default=False,
        help="Keep sites with zero difference (default: filter them).",
    )
    qc_group.add_argument(
        "--disable-preprocess-two-sample",
        action="store_true",
        default=False,
        help="Disable preprocessing for two-sample tests.",
    )
    # Tiered Resource Allocation Defaults
    resource_group.add_argument(
        "--threads-per-job",
        type=int,
        default=16,
        help="Default threads for a single parallelizable job.",
    )
    resource_group.add_argument(
        "--low-mem",
        type=int,
        default=10000,
        help="Memory (MB) for low-memory jobs (e.g., profile, significance_test).",
    )
    resource_group.add_argument(
        "--high-mem",
        type=int,
        default=100000,
        help="Memory (MB) for high-memory jobs (e.g., quality_control, analyze_alleles).",
    )
    resource_group.add_argument(
        "--low-time",
        type=str,
        default="08:00:00",
        help="Runtime for low-runtime jobs (e.g., profile, significance_test).",
    )
    resource_group.add_argument(
        "--high-time",
        type=str,
        default="24:00:00",
        help="Runtime for high-runtime jobs (e.g., general tasks).",
    )

    args = parser.parse_args()

    # --- Argument Validation ---
    data_type_from_args = args.data_type
    # take data_type from config if not provided in args
    if not data_type_from_args and args.config:
        with open(args.config, "r") as f:
            data_type_from_args = yaml.safe_load(f).get("analysis", {}).get("data_type")

    if args.timepoints:
        if data_type_from_args == "longitudinal":
            if len(args.timepoints) != 2:
                parser.error(
                    "--timepoints requires exactly 2 values in longitudinal mode."
                )
            if not args.focus:
                parser.error("--focus is required in longitudinal mode.")
            if args.focus not in args.timepoints:
                parser.error(
                    f"--focus value '{args.focus}' must be one of the provided timepoints: {args.timepoints}"
                )
        elif data_type_from_args == "single":
            if len(args.timepoints) != 1:
                parser.error("--timepoints requires exactly 1 value in single mode.")
            if args.focus:
                parser.error("--focus cannot be used in single mode.")

    runtime_config_path = (
        f"runtime_config_{datetime.now().strftime('%Y_%m_%d-%H_%M_%S')}.yml"
    )

    config_data = create_config(args)

    # --- Generate Config and Exit (if requested) ---
    config_data = create_config(args)
    if args.generate_config:
        with open(args.generate_config, "w") as f:
            yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
        logging.info(
            f"Successfully generated configuration file at: {args.generate_config}"
        )
        sys.exit(0)

    # --- Normal Workflow Execution ---
    with open(runtime_config_path, "w") as f:
        yaml.dump(config_data, f, default_flow_style=False, sort_keys=False)
    logging.info(f"Created runtime config file for this run: {runtime_config_path}")

    runner_resource = resources.files("alleleflux").joinpath("runner.sh")

    with resources.as_file(runner_resource) as runner_script_path:
        runner_script_path.chmod(0o755)

        command = [str(runner_script_path), "--config", runtime_config_path]

        if args.cores:
            command.extend(["--cores", str(args.cores)])
        if args.profile:
            command.extend(["--profile", args.profile])
        if args.dry_run:
            command.append("--dry-run")
        if args.step1_only:
            command.append("--step1-only")
        if args.step2_only:
            command.append("--step2-only")

        env = os.environ.copy()
        env["ALLELEFLUX_VERSION"] = metadata.version("AlleleFlux")

        process = subprocess.Popen(command, env=env, preexec_fn=os.setsid)

        try:
            exit_code = process.wait()
            sys.exit(exit_code)
        except KeyboardInterrupt:
            logging.warning(
                "\nKeyboard interrupt received. Sending SIGINT to child process group..."
            )
            os.killpg(os.getpgid(process.pid), signal.SIGINT)
            process.wait()
            sys.exit(130)


if __name__ == "__main__":
    main()
