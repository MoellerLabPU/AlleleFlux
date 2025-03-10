#!/bin/bash

# Define snakefile paths for the two workflows.
SNAKEFILE1="smk_workflow/snakefile.smk"
SNAKEFILE2="smk_workflow/snakefile2.smk"

PROFILE="smk_workflow/della_profile/"

# Run Workflow 1: Preprocessing
echo "Running Workflow 1 (Preprocessing) with profile '$PROFILE'"
snakemake --snakefile "$SNAKEFILE1" --profile "$PROFILE"
if [ $? -ne 0 ]; then
    echo "Workflow 1 failed! Aborting."
    exit 1
fi


# Run Workflow 2: Allele Analysis
echo "Running Workflow 2 (Allele Analysis) with profile '$PROFILE'"
snakemake --snakefile "$SNAKEFILE2" --profile "$PROFILE"
if [ $? -ne 0 ]; then
    echo "Workflow 2 failed! Aborting."
    exit 1
fi

echo "Both workflows completed successfully!"