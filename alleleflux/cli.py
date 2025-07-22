#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

import click
import pkg_resources
import questionary
import yaml

from . import __version__


# Helper function to process boolean inputs with defaults
def process_boolean_input(input_value, default):
    """Process a text input as a boolean value."""
    if not input_value:
        return default

    input_value = input_value.lower().strip()
    if input_value in ["true", "t", "yes", "y", "1"]:
        return True
    elif input_value in ["false", "f", "no", "n", "0"]:
        return False
    else:
        click.echo(f"Invalid value '{input_value}'. Using default: {default}")
        return default


# Helper to get the path to the template file
def get_template_path():
    return pkg_resources.resource_filename("alleleflux", "config.template.yml")


@click.group()
def cli():
    """AlleleFlux: A tool for analyzing allele frequency changes in metagenomic data."""
    pass


@cli.command()
@click.option(
    "--template",
    "use_template",
    is_flag=True,
    help="Print a template config to stdout.",
)
@click.option(
    "--output",
    default="alleleflux_config.yml",
    help="Output config file path.",
    show_default=True,
)
def init_config(use_template, output):
    """Creates a new configuration file for an AlleleFlux analysis."""
    template_path = get_template_path()
    with open(template_path, "r") as f:
        config_template = yaml.safe_load(f)

    if use_template:
        with open(template_path, "r") as f:
            click.echo(f.read())
        return

    click.echo("Welcome to AlleleFlux! Let's create your configuration file.\n")

    # --- Interactive Prompts ---
    try:
        # Required input files
        click.echo("📁 Required Input Files:")
        bam_dir = questionary.text(
            "What is the path to your directory of BAM files?",
            validate=lambda text: (
                True
                if Path(text).is_dir()
                else "Path does not exist or is not a directory."
            ),
        ).ask()
        if bam_dir is None:
            raise KeyboardInterrupt

        fasta_path = questionary.text(
            "What is the path to your reference FASTA file?",
            validate=lambda text: (
                True if Path(text).is_file() else "File does not exist."
            ),
        ).ask()
        if fasta_path is None:
            raise KeyboardInterrupt

        metadata_path = questionary.text(
            "What is the path to your metadata file?",
            validate=lambda text: (
                True if Path(text).is_file() else "File does not exist."
            ),
        ).ask()
        if metadata_path is None:
            raise KeyboardInterrupt

        prodigal_path = questionary.text(
            "Path to Prodigal genes file (.fna):",
            validate=lambda text: (
                True if Path(text).is_file() else "File does not exist."
            ),
        ).ask()
        if prodigal_path is None:
            raise KeyboardInterrupt

        gtdb_path = questionary.text(
            "Path to GTDB taxonomy file:",
            validate=lambda text: (
                True if Path(text).is_file() else "File does not exist."
            ),
        ).ask()
        if gtdb_path is None:
            raise KeyboardInterrupt

        mag_mapping_path = questionary.text(
            "Path to MAG mapping file:",
            validate=lambda text: (
                True if Path(text).is_file() else "File does not exist."
            ),
        ).ask()
        if mag_mapping_path is None:
            raise KeyboardInterrupt

        # Output directory
        click.echo("\n📂 Output Configuration:")
        output_dir = questionary.text(
            "Where should the output files be saved?", default="./alleleflux_output"
        ).ask()
        if output_dir is None:
            raise KeyboardInterrupt

        # Analysis type
        click.echo("\n🔬 Analysis Configuration (press Enter to use defaults):")

        # Default values for analysis configuration
        default_data_type = "single"
        default_use_lmm = True
        default_use_significance_tests = True
        default_use_cmh = False

        # Ask for data type with explanation
        click.echo("Analysis Type Options:")
        click.echo("  • single: For analyzing allele frequencies at a single timepoint")
        click.echo(
            "  • longitudinal: For analyzing changes in allele frequencies across multiple timepoints"
        )

        data_type_input = questionary.text(
            "Type of analysis [single]:", default="single"
        ).ask()
        if data_type_input is None:
            raise KeyboardInterrupt

        # Process and validate data type
        data_type = (
            data_type_input.lower().strip() if data_type_input else default_data_type
        )
        if data_type not in ["single", "longitudinal"]:
            click.echo(
                f"Invalid data type '{data_type}'. Using default: {default_data_type}"
            )
            data_type = default_data_type

        # Update CMH default based on data type
        if data_type == "longitudinal":
            default_use_cmh = True

        # Ask for remaining analysis options with clear descriptions
        click.echo("\nStatistical Analysis Options:")

        use_lmm_input = questionary.text(
            "Use Linear Mixed Models analysis? (true/false) [true]:", default=""
        ).ask()
        if use_lmm_input is None:
            raise KeyboardInterrupt
        use_lmm = process_boolean_input(use_lmm_input, default_use_lmm)

        use_significance_tests_input = questionary.text(
            "Use significance tests? (true/false) [true]:", default=""
        ).ask()
        if use_significance_tests_input is None:
            raise KeyboardInterrupt
        use_significance_tests = process_boolean_input(
            use_significance_tests_input, default_use_significance_tests
        )

        use_cmh_input = questionary.text(
            f"Use Cochran-Mantel-Haenszel test? (true/false) [{str(default_use_cmh).lower()}]:",
            default="",
        ).ask()
        if use_cmh_input is None:
            raise KeyboardInterrupt
        use_cmh = process_boolean_input(use_cmh_input, default_use_cmh)

        # Timepoint configuration
        click.echo("\n⏱️  Timepoint Configuration:")

        if data_type == "single":
            click.echo(
                "For single timepoint analysis, please specify the timepoint name(s)."
            )
            click.echo(
                "You can add multiple timepoints by entering them one at a time."
            )

            timepoints = []
            while True:
                timepoint = questionary.text(
                    "Enter a timepoint name (or leave blank to finish):",
                ).ask()
                if timepoint is None:
                    raise KeyboardInterrupt

                if not timepoint:
                    if not timepoints:
                        click.echo(
                            "At least one timepoint is required. Using default 'end'"
                        )
                        timepoints.append("end")
                    break

                timepoints.append(timepoint)

            # Create timepoints_combinations for single timepoint analysis
            timepoints_combinations = []
            for tp in timepoints:
                timepoints_combinations.append({"timepoint": [tp]})

        else:  # longitudinal
            click.echo(
                "For longitudinal analysis, please specify timepoint pairs and focus timepoint."
            )

            timepoints_combinations = []
            while True:
                # First timepoint
                timepoint1 = questionary.text(
                    "Enter first timepoint name (or leave blank to finish):",
                ).ask()
                if timepoint1 is None:
                    raise KeyboardInterrupt

                if not timepoint1:
                    if not timepoints_combinations:
                        click.echo(
                            "At least one timepoint pair is required. Using default 'pre' and 'post'"
                        )
                        timepoints_combinations.append(
                            {"timepoint": ["pre", "post"], "focus": "post"}
                        )
                    break

                # Second timepoint
                timepoint2 = questionary.text(
                    "Enter second timepoint name:",
                    validate=lambda text: (
                        True
                        if text and text != timepoint1
                        else "Please enter a different timepoint name."
                    ),
                ).ask()
                if timepoint2 is None:
                    raise KeyboardInterrupt

                # Focus timepoint
                focus = questionary.select(
                    "Select the focus timepoint:",
                    choices=[
                        {"name": timepoint1, "value": timepoint1},
                        {"name": timepoint2, "value": timepoint2},
                    ],
                ).ask()
                if focus is None:
                    raise KeyboardInterrupt

                timepoints_combinations.append(
                    {"timepoint": [timepoint1, timepoint2], "focus": focus}
                )

        # Group configuration
        click.echo("\n👥 Group Configuration:")
        click.echo("Please specify groups to compare (e.g., 'treatment', 'control').")

        groups = []
        while True:
            group = questionary.text(
                "Enter a group name (or leave blank to finish):",
            ).ask()
            if group is None:
                raise KeyboardInterrupt

            if not group:
                if len(groups) < 2:
                    click.echo(
                        "At least two groups are required. Using default 'treatment' and 'control'"
                    )
                    groups = ["treatment", "control"]
                break

            groups.append(group)

            if len(groups) >= 2:
                add_more = questionary.confirm(
                    "Add another group?", default=False
                ).ask()
                if add_more is None:
                    raise KeyboardInterrupt
                if not add_more:
                    break

        # Quality control parameters
        click.echo("\n⚙️  Quality Control Parameters (press Enter to use defaults):")

        min_sample_num_input = questionary.text(
            "Minimum number of samples required for analysis [4]:", default=""
        ).ask()
        if min_sample_num_input is None:
            raise KeyboardInterrupt

        # Process and validate min_sample_num
        if min_sample_num_input:
            if min_sample_num_input.isdigit() and int(min_sample_num_input) > 0:
                min_sample_num = int(min_sample_num_input)
            else:
                click.echo(f"Invalid value '{min_sample_num_input}'. Using default: 4")
                min_sample_num = 4
        else:
            min_sample_num = 4

        breadth_threshold_input = questionary.text(
            "Minimum coverage breadth threshold (0-1) [0.1]:", default=""
        ).ask()
        if breadth_threshold_input is None:
            raise KeyboardInterrupt

        # Process and validate breadth_threshold
        if breadth_threshold_input:
            try:
                bt_value = float(breadth_threshold_input)
                if 0 <= bt_value <= 1:
                    breadth_threshold = bt_value
                else:
                    click.echo(
                        f"Invalid value '{breadth_threshold_input}'. Using default: 0.1"
                    )
                    breadth_threshold = 0.1
            except ValueError:
                click.echo(
                    f"Invalid value '{breadth_threshold_input}'. Using default: 0.1"
                )
                breadth_threshold = 0.1
        else:
            breadth_threshold = 0.1

    except KeyboardInterrupt:
        click.echo("\nConfiguration cancelled.")
        return

    # --- Populate the config ---
    config_template["input"]["bam_dir"] = bam_dir
    config_template["input"]["fasta_path"] = fasta_path
    config_template["input"]["metadata_path"] = metadata_path
    config_template["input"]["prodigal_path"] = prodigal_path
    config_template["input"]["gtdb_path"] = gtdb_path
    config_template["input"]["mag_mapping_path"] = mag_mapping_path

    config_template["output"]["root_dir"] = output_dir

    # Analysis configuration
    config_template["analysis"]["data_type"] = data_type
    config_template["analysis"]["use_lmm"] = use_lmm
    config_template["analysis"]["use_significance_tests"] = use_significance_tests
    config_template["analysis"]["use_cmh"] = use_cmh

    # Quality control
    config_template["quality_control"]["min_sample_num"] = int(min_sample_num)
    config_template["quality_control"]["breadth_threshold"] = float(breadth_threshold)

    # Set timepoints and groups
    config_template["analysis"]["timepoints_combinations"] = timepoints_combinations
    config_template["analysis"]["groups_combinations"] = [groups]

    # --- Save the new config file ---
    with open(output, "w") as f:
        yaml.dump(config_template, f, sort_keys=False, default_flow_style=False)

    click.echo(f"\n✅ Configuration file '{output}' was created successfully!")
    click.echo("\n🚀 To start your analysis, run the following command:")
    click.echo(
        f"   alleleflux run --config {output}"
    )  # Assuming a 'run' command will exist


if __name__ == "__main__":
    cli()
