"""
dN/dS Analysis Rules
Rules for performing dN/dS analysis from timepoints using dnds_from_timepoints.py
"""

def _get_dnds_test_type_str(test_type, for_eligibility=False):
    """Helper function to get the standardized test type string for dN/dS analysis."""
    if test_type in ["two_sample_unpaired_tTest", "two_sample_unpaired_MannWhitney"]:
        return "two_sample_unpaired"
    elif test_type in ["two_sample_paired_tTest", "two_sample_paired_Wilcoxon"]:
        return "two_sample_paired"
    elif test_type in ["single_sample_tTest", "single_sample_Wilcoxon"]:
        return "single_sample"
    elif test_type in ["lmm", "lmm_across_time", "cmh", "cmh_across_time"]:
        return test_type
    else:
        raise ValueError(f"Unsupported DN_DS_TEST_TYPE: {test_type}")

def get_significant_sites_df_path():
    """Get the path to the significant sites file based on the test type."""
    test_type_str = _get_dnds_test_type_str(DN_DS_TEST_TYPE)
    return os.path.join(
        OUTDIR,
        "p_value_summary",
        "{timepoints}-{groups}",
        f"p_value_summary_{test_type_str}_{{timepoints}}.tsv"
    )

rule dnds_from_timepoints:
    """
    Run dN/dS analysis for all eligible MAGs for a given subject pair.
    This rule executes dnds_from_timepoints.py once per subject, passing all
    eligible MAGs, and outputs to a directory.
    """
    input:
        # P-value summary file for the specified test type
        significant_sites=get_significant_sites_df_path(),
        # Profile directory containing all sample profiles
        profile_dir=os.path.join(OUTDIR, "profiles"),
        # Prodigal gene predictions
        prodigal_fasta=config["input"]["prodigal_path"]
    output:
        directory(
            os.path.join(
                OUTDIR,
                "dnds_analysis",
                "{timepoints}-{groups}",
                "{subject_id}"
            )
        )
    params:
        p_value_column=config["dnds"]["p_value_column"],
        p_value_threshold=config["statistics"]["p_value_threshold"],  # Use from statistics section
        dn_ds_test_type=DN_DS_TEST_TYPE,
    resources:
        time=config["resources"]["time"]["general"],
    threads: config["resources"]["cpus"]["threads_per_job"]
    run:
        import os
        import pandas as pd
        
        # Determine the eligibility test type string based on the specific test
        eligibility_test_type = _get_dnds_test_type_str(params.dn_ds_test_type, for_eligibility=True)

        # Get all eligible MAGs for this specific test type
        eligibility_file = os.path.join(
            OUTDIR, f"eligibility_table_{wildcards.timepoints}-{wildcards.groups}.tsv"
        )
        df = pd.read_csv(eligibility_file, sep="\t")

        if eligibility_test_type in ["two_sample_unpaired", "lmm"]:
            eligible_mags = df.loc[df["unpaired_test_eligible"] == True, "MAG_ID"].tolist()
        elif eligibility_test_type in ["two_sample_paired", "cmh"]:
            eligible_mags = df.loc[df["paired_test_eligible"] == True, "MAG_ID"].tolist()
        elif eligibility_test_type in ["single_sample", "lmm_across_time", "cmh_across_time"]:
            single_cols = [col for col in df.columns if col.startswith("single_sample_eligible_")]
            if not single_cols:
                eligible_mags = []
            else:
                eligible_mags = df.loc[df[single_cols].any(axis=1), "MAG_ID"].unique().tolist()
        else:
            raise ValueError(f"Unknown eligibility type for dN/dS: {eligibility_test_type}")

        if not eligible_mags:
            # If no MAGs are eligible, create an empty output directory and touch a sentinel file
            shell(f"mkdir -p {output} && touch {output}/no_eligible_mags")
            print(f"No eligible MAGs for {wildcards.timepoints}-{wildcards.groups}. Created empty directory.")
            # Stop further execution of the rule
            return

        mag_ids_str = " ".join(eligible_mags)
        
        # Parse sample pairs for this timepoint-group combination
        sample_pairs = parse_metadata_for_timepoint_pairs(wildcards.timepoints, wildcards.groups)
        subject_pair = None
        for subject_id, ancestral_sample_id, derived_sample_id in sample_pairs:
            if str(subject_id) == wildcards.subject_id:
                subject_pair = (ancestral_sample_id, derived_sample_id)
                break
        
        if subject_pair is None:
            raise ValueError(f"No sample pair found for subject {wildcards.subject_id}")
        
        ancestral_sample_id, derived_sample_id = subject_pair

        # Conditionally add the --group-analyzed flag
        group_analyzed_flag = ""
        if eligibility_test_type in ["single_sample", "lmm_across_time", "cmh_across_time"]:
            # Ensure the group_analyzed key is present in the config
            if "group_analyzed" not in config["dnds"]:
                raise ValueError(
                    f"The test type '{params.dn_ds_test_type}' requires 'group_analyzed' to be set in the dnds section of the config."
                )
            group = config["dnds"]["group_analyzed"]
            group_analyzed_flag = f"--group-analyzed {group}"

        shell(
            """
            alleleflux-dnds-from-timepoints \
                --significant_sites {input.significant_sites} \
                --mag_ids {mag_ids_str} \
                --p_value_column {params.p_value_column} \
                --p_value_threshold {params.p_value_threshold} \
                --test-type {params.dn_ds_test_type} \
                {group_analyzed_flag} \
                --ancestral_sample_id {ancestral_sample_id} \
                --derived_sample_id {derived_sample_id} \
                --profile_dir {input.profile_dir} \
                --prodigal_fasta {input.prodigal_fasta} \
                --outdir {output} \
                --prefix {wildcards.subject_id} \
                --cpus {threads} 
            """
        )
