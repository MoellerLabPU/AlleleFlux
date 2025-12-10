"""Within-group significance testing rules.

This module contains rules for statistical tests analyzing allele frequency
changes within individual experimental groups across timepoints:
- Single-sample tests (one-sample t-test, Wilcoxon signed-rank)
- LMM across-time analysis within groups
- CMH across-time analysis within groups
- Preprocessing rules for filtering zero-count positions
"""


rule preprocess_within_groups:
    input:
        mean_allele_changes=os.path.join(
            OUTDIR,
            "allele_analysis",
            "allele_analysis_{timepoints}-{groups}",
            "{mag}_allele_frequency_changes_mean.tsv.gz",
        ),
    output:
        outPath=get_preprocessed_within_groups_path(),
        statusPath=os.path.join(
            OUTDIR,
            "significance_tests",
            "preprocessed_within_groups_{timepoints}-{groups}",
            "{mag}_{group}_preprocessing_status.json",
        ),
    params:
        max_zero_flag=(
            ("--max_zero_count " + str(config["statistics"]["max_zero_count"]))
            if config["statistics"].get("max_zero_count", None) is not None
            else ""
        ),
        outDir=os.path.join(
            OUTDIR, "significance_tests", "preprocessed_within_groups_{timepoints}-{groups}"
        ),
        min_positions=config["statistics"].get("min_positions_after_preprocess", 1),
        min_sample_num=config["quality_control"]["min_sample_num"],
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-preprocess-within-group \
            --df_fPath {input.mean_allele_changes} \
            --mag_id {wildcards.mag} \
            --group {wildcards.group} \
            --output_dir {params.outDir} \
            --min_positions {params.min_positions} \
            --min_sample_num {params.min_sample_num} \
            {params.max_zero_flag}
            """

rule single_sample:
    input:
        mean_allele_changes=(
            get_preprocessed_within_groups_path()
            if config["statistics"].get("preprocess_within_groups", False)
            else os.path.join(
                OUTDIR,
                "allele_analysis",
                "allele_analysis_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_mean.tsv.gz",
            )
        ),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "single_sample_{timepoints}-{groups}",
            "{mag}_single_sample_{group}.tsv.gz",
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "single_sample_{timepoints}-{groups}"
        ),
    threads: config["resources"]["cpus"]["threads_per_job"]
    resources:
        time=config["resources"]["time"]["significance_test"],
    shell:
        """
        alleleflux-single-sample \
            --df_fPath {input.mean_allele_changes} \
            --min_sample_num {params.min_sample_num} \
            --mag_id {wildcards.mag} \
            --group {wildcards.group} \
            --cpus {threads} \
            --output_dir {params.outDir} 
            """


rule lmm_analysis_across_time:
    input:
        input_df=(
            os.path.join(
                OUTDIR,
                "allele_analysis",
                "allele_analysis_{timepoints}-{groups}",
                "{mag}_allele_frequency_changes_no_zero-diff.tsv.gz",
            )
            if not config["quality_control"].get("disable_zero_diff_filtering", False)
            else get_longitudinal_input_path()
        ),
        preprocessed_df=(
            get_preprocessed_within_groups_path()
            if config["statistics"].get("preprocess_within_groups", False)
            else []
        ),
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "lmm_across_time_{timepoints}-{groups}",
            "{mag}_lmm_across_time_{group}.tsv.gz"
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "lmm_across_time_{timepoints}-{groups}"
        ),
        input_time_format=(
            "suffix"
            if not config["quality_control"].get("disable_zero_diff_filtering", False)
            else "column"
        ),
        preprocessed_flag=(
            "--preprocessed_df"
            if config["statistics"].get("preprocess_within_groups", False)
            else ""
        )
    threads: config["resources"]["cpus"]["threads_per_job"]
    resources:
        time=config["resources"]["time"]["significance_test"],
        mem_mb=config["resources"]["memory"]["significance_test"] # LMM can be memory intensive
    shell:
        """
        alleleflux-lmm \
            --input_df {input.input_df} \
            --mag_id {wildcards.mag} \
            --output_dir {params.outDir} \
            --data_type across_time \
            --group_to_analyze {wildcards.group} \
            --input_time_format {params.input_time_format} \
            --cpus {threads} \
            --min_sample_num {params.min_sample_num} \
            {params.preprocessed_flag} {input.preprocessed_df}
        """


rule cmh_test_across_time:
    input:
        input_df=get_longitudinal_input_path(),
        preprocessed_df=(
            get_preprocessed_within_groups_path()
            if config["statistics"].get("preprocess_within_groups", False)
            else []
        )
    output:
        os.path.join(
            OUTDIR,
            "significance_tests",
            "cmh_across_time_{timepoints}-{groups}",
            "{mag}_cmh_across_time_{group}.tsv.gz"
        ),
    params:
        min_sample_num=config["quality_control"]["min_sample_num"],
        outDir=os.path.join(
            OUTDIR, "significance_tests", "cmh_across_time_{timepoints}-{groups}"
        ),
        preprocessed_flag=(
            "--preprocessed_df"
            if config["statistics"].get("preprocess_within_groups", False)
            else ""
        )
    threads: config["resources"]["cpus"]["threads_per_job"]
    resources:
        time=config["resources"]["time"]["significance_test"],
        mem_mb=config["resources"]["memory"]["significance_test"],
    shell:
        """
        alleleflux-cmh \
            --input_df {input.input_df} \
            --mag_id {wildcards.mag} \
            --output_dir {params.outDir} \
            --data_type across_time \
            --group {wildcards.group} \
            --cpus {threads} \
            --min_sample_num {params.min_sample_num} \
            {params.preprocessed_flag} {input.preprocessed_df}
        """
