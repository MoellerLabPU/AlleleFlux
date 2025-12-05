# Snakefile
# Description: A Snakemake workflow to clean and combine metadata tables.
METADATA_OUT = os.path.join(config["output_dir"], "combined_metadata.tsv")

# Create a dictionary to map run IDs to their full parameter sets
# e.g., {"run1": {"name": "run1", "path": ...}, "run2": ...}
RUN_ID_PARAMS = {run_id['name']: run_id for run_id in config["metadata_params"]}

# Get the list of run IDs from the config
# e.g., ["run1", "run2"]
RUN_IDS = list(RUN_ID_PARAMS.keys())


localrules: 
    clean_metadata,
    aggregate_metadata,

# === RULE : 'clean_metadata' (The "Fan-out" Step) ===
# This rule takes one raw metadata file, cleans it using your
# script, and produces one corresponding cleaned file.
# This rule will be run once for every item in the 'RUN_IDS' list.
rule clean_metadata:
    input:
        # Input file is now dynamically fetched from the config
        tsv_in=lambda wildcards: RUN_ID_PARAMS[wildcards.run_id]["path"]
    output:
        # Output file, e.g., "results/cleaned/run1.tsv"
        # 'temp' marks it for deletion after the final file is made.
        tsv_out=temp(os.path.join(config["output_dir"], "cleaned_metadata", "{run_id}.tsv"))
    params:
        # Pass the entire parameter set for the specific run
        run_params=lambda wildcards: RUN_ID_PARAMS[wildcards.run_id]
    run:
        # Base command
        cmd = [
            "alleleflux-prepare-metadata",
            f"--metadata-in {input.tsv_in}",
            f"--metadata-out {output.tsv_out}",
            f"--base-profile-dir '{params.run_params['base_profile_dir']}'"
        ]

        # Add columns with defaults
        cmd.append(f"--subject-col '{params.run_params.get('subject_col', 'subjectID')}'")
        cmd.append(f"--sample-col '{params.run_params.get('sample_col', 'sample_id')}'")
        cmd.append(f"--group-col '{params.run_params.get('group_col', 'group')}'")
        cmd.append(f"--time-col '{params.run_params.get('time_col', 'time')}'")
        cmd.append(f"--day-col '{params.run_params.get('day_col', 'day')}'")
        cmd.append(f"--replicate-col '{params.run_params.get('replicate_col', 'replicate')}'")

        # Execute the command
        shell(" ".join(cmd))

# === RULE : 'aggregate_metadata' (The "Gather" Step) ===
# This rule runs *after* all 'clean_metadata' jobs are finished.
# It "gathers" all the individual cleaned files and concatenates
# them into the one final table.
rule aggregate_metadata:
    input:
        # Gather all cleaned files, e.g.,
        # ["results/cleaned/run1.tsv", "results/cleaned/run2.tsv", ...]
        cleaned_files=expand(os.path.join(config["output_dir"], "cleaned_metadata", "{run_id}.tsv"), run_id=RUN_IDS)
    output:
        # The final target file, from the top-level 'metadata_out' key
        final_file=METADATA_OUT
    run:
        import pandas as pd
        from snakemake.logging import logger

        logger.info(f"Aggregating {len(input.cleaned_files)} files...")
        
        # Read all intermediate files into a list of DataFrames
        all_dfs = [pd.read_csv(file, sep='\t') for file in input.cleaned_files]
        
        # Combine all DataFrames into one
        combined_df = pd.concat(all_dfs, ignore_index=True)
        
        # Check for duplicates and raise error if any exist
        if combined_df.duplicated().any():
            raise ValueError("Duplicate rows found in the combined metadata. Please resolve duplicates before proceeding.")
        
        # Save the final, combined table
        combined_df.to_csv(str(output.final_file), sep="\t", index=False)
        
        logger.info(f"Success! Saved {len(combined_df)} total rows to {output.final_file}")