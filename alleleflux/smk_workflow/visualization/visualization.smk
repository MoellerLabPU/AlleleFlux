import os
import pandas as pd

# Load configuration
configfile: "visualization_config.yaml"

# Include sub-workflows
include: "smks/metadata_prep.smk"
include: "smks/analysis.smk"
include: "smks/plotting.smk"

# Main rule to drive the workflow
rule all:
    input:
        METADATA_OUT,
        get_tracking_targets,
        get_plotting_targets
