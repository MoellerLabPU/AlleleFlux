"""dN/dS (nonsynonymous to synonymous substitution ratio) analysis rules.

This module performs dN/dS analysis for longitudinal data by comparing
ancestral and derived samples from the same subject. It identifies
evolutionary patterns by calculating substitution rates at significant
allele frequency change positions.
"""


def _get_dnds_test_type_str(test_type, for_eligibility=False):
    """Helper function to get the standardized test type string for dN/dS analysis."""
    if test_type in ["two_sample_unpaired_tTest", "two_sample_unpaired_MannWhitney", "two_sample_unpaired_tTest_abs", "two_sample_unpaired_MannWhitney_abs"]:
        return "two_sample_unpaired"
    elif test_type in ["two_sample_paired_tTest", "two_sample_paired_Wilcoxon", "two_sample_paired_tTest_abs", "two_sample_paired_Wilcoxon_abs"]:
        return "two_sample_paired"
    elif test_type in ["single_sample_tTest", "single_sample_Wilcoxon"]:
        return "single_sample"
    elif test_type in ["lmm", "lmm_abs", "lmm_across_time", "cmh", "cmh_across_time"]:
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
        # Uses PROFILES_DIR which points to either existing profiles (if profiles_path
        # is specified in config) or newly generated profiles in OUTDIR/profiles
        profile_dir=PROFILES_DIR,
        # Prodigal gene predictions
        prodigal_fasta=config["input"]["prodigal_path"],
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
        log_level=config.get("log_level", "INFO"),
    resources:
        time=get_time("dnds_from_timepoints"),
        mem_mb=get_mem_mb("dnds_from_timepoints"),
    threads: get_threads("dnds_from_timepoints")
    run:
        import os
        import subprocess
        from pathlib import Path
        import pandas as pd
        from snakemake.logging import logger
        # from alleleflux.scripts.utilities.logging_config import setup_logging
        # setup_logging()
        # logger = logging.getLogger(__name__)

        outdir = Path(str(output[0]))
        
        valid_p_value_columns = ["min_p_value", "q_value"]
        if params.p_value_column not in valid_p_value_columns:
            raise ValueError(
                f"Invalid p_value_column '{params.p_value_column}'. "
                f"Must be one of: {', '.join(valid_p_value_columns)}"
            )
        
        # Load the p-value summary to find MAGs with at least one significant site
        sig_sites_df = pd.read_csv(input.significant_sites, sep="\t")
        if params.p_value_column not in sig_sites_df.columns:
            raise ValueError(
                f"P-value column '{params.p_value_column}' not found in {input.significant_sites}. "
                f"Available columns: {', '.join(sig_sites_df.columns)}"
            )

        # If group_analyzed is specified, filter the significant sites dataframe by that group first
        group_analyzed_flag = ""
        if params.dn_ds_test_type in ["single_sample_tTest", "single_sample_Wilcoxon", "lmm_across_time", "cmh_across_time"]:
            group = config["dnds"].get("group_analyzed")
            if not group:
                raise ValueError(
                    f"For test type '{params.dn_ds_test_type}', 'group_analyzed' must be specified in the 'dnds' section of the config."
                )
            if "group_analyzed" not in sig_sites_df.columns:
                raise ValueError(
                    f"A 'group_analyzed' column is expected in {input.significant_sites} for test type "
                    f"'{params.dn_ds_test_type}', but is not present."
                )
            # Conditionally add the --group-analyzed flag
            group_analyzed_flag = f"--group-analyzed {group}"
            
            sig_sites_df = sig_sites_df[sig_sites_df["group_analyzed"] == group]
            logger.info(f"Filtered significant sites for group '{group}'.")

        sig_sites_df = sig_sites_df[sig_sites_df[params.p_value_column] <= params.p_value_threshold]

        test_type = params.dn_ds_test_type
        if params.dn_ds_test_type in ["lmm", "lmm_across_time"]:
            test_type = "LMM"
        elif params.dn_ds_test_type in ["cmh", "cmh_across_time"]:
            test_type = "CMH"

        sig_sites_df = sig_sites_df[sig_sites_df["test_type"] == test_type]
        eligible_mags = sig_sites_df["mag_id"].unique().tolist()
        
        if eligible_mags:
            logger.info(
                f"Found {len(eligible_mags)} MAG(s) that are both eligible and have significant sites "
                f"({params.p_value_column} <= {params.p_value_threshold})."
            )


        if not eligible_mags:
            # If no MAGs are eligible, create an empty output directory and touch a sentinel file
            outdir.mkdir(parents=True, exist_ok=True)
            (outdir / "no_eligible_mags").touch()
            logger.info(f"No eligible MAGs for {wildcards.timepoints}-{wildcards.groups}. Created empty directory.")
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
        
        # Log which samples are ancestral and derived
        logger.info(
            f"dN/dS analysis for subject {wildcards.subject_id}: "
            f"Ancestral (Time 1) = {ancestral_sample_id}, "
            f"Derived (Time 2) = {derived_sample_id}"
        )
        
        cmd = f"""
        alleleflux-dnds-from-timepoints \\
            --significant_sites {input.significant_sites} \\
            --mag_ids {mag_ids_str} \\
            --p_value_column {params.p_value_column} \\
            --p_value_threshold {params.p_value_threshold} \\
            --test-type {test_type} \\
            {group_analyzed_flag} \\
            --ancestral_sample_id {ancestral_sample_id} \\
            --derived_sample_id {derived_sample_id} \\
            --profile_dir {input.profile_dir} \\
            --prodigal_fasta {input.prodigal_fasta} \\
            --outdir {output} \\
            --prefix {wildcards.subject_id} \\
            --cpus {threads} \\
            --log-level {params.log_level}
        """

        # Bypass Snakemake's shell() logging which is causing TypeError
        logger.info(f"Executing command: {cmd}")
        
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed with exit code {e.returncode}")
            raise e
