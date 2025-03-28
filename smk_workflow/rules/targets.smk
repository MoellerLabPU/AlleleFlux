def get_allele_analysis_targets():
    """
    Dynamically generate targets for allele analysis based on eligible MAGs.
    """
    targets = []
    for timepoints in timepoints_labels:
        for groups in groups_labels:
            # Get eligible MAGs for this timepoint-group combination
            eligible_mags = get_mags_by_eligibility(timepoints, groups, eligibility_type="all")
            # Add targets for each eligible MAG
            for mag in eligible_mags:
                # Different output files based on data_type and filtering options
                if DATA_TYPE == "single":
                    if not config["quality_control"].get("disable_zero_diff_filtering", False):
                        # When single data type and filtering is not disabled
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "allele_analysis",
                                f"allele_analysis_{timepoints}-{groups}",
                                f"{mag}_allele_frequency_no_constant.tsv.gz",
                            )
                        )
                    else:
                        # When single data type and filtering is disabled
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "allele_analysis",
                                f"allele_analysis_{timepoints}-{groups}",
                                f"{mag}_allele_frequency_single.tsv.gz",
                            )
                        )
                else:  # longitudinal
                    targets.append(
                        os.path.join(
                            OUTDIR,
                            "allele_analysis",
                            f"allele_analysis_{timepoints}-{groups}",
                            f"{mag}_allele_frequency_changes_mean.tsv.gz",
                        )
                    )
    return targets


def get_taxa_scores_targets():
    targets = []
    tax_levels = ["phylum", "class", "order", "family", "genus", "species"]
    
    # For two-sample tests, group_str is empty.
    if config["analysis"].get("use_significance_tests", True):
        for tp in timepoints_labels:
            for gr in groups_labels:
                for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                    
                    # Only generate targets if there are eligible MAGs for this test type
                    mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                    if mags:  # Only proceed if there are eligible MAGs
                        group_str = ""  # no group marker for two-sample tests
                        for taxon in tax_levels:
                            targets.append(
                                os.path.join(
                                    OUTDIR,
                                    "scores",
                                    "processed",
                                    "combined",
                                    f"scores_{test_type}-{tp}-{gr}{group_str}-{taxon}.tsv",
                                )
                            )
                    # print(f"Targets for unpaired are: {mags}")

        
        # For single-sample tests, group_str is "_" plus the sample group.
        # Only include if data_type is longitudinal
        if DATA_TYPE == "longitudinal":
            for tp in timepoints_labels:
                for gr in groups_labels:
                    sample_entries = get_single_sample_entries(tp, gr)
                    if sample_entries:  # Only proceed if there are eligible MAGs
                        unique_groups = sorted(set([grp for mag, grp in sample_entries]))
                        for grp in unique_groups:
                            group_str = f"_{grp}"
                            for taxon in tax_levels:
                                targets.append(
                                    os.path.join(
                                        OUTDIR,
                                        "scores",
                                        "processed",
                                        "combined",
                                        f"scores_single_sample-{tp}-{gr}{group_str}-{taxon}.tsv",
                                    )
                                )
    
    # Add LMM taxa targets if LMM is enabled
    if config["analysis"].get("use_lmm", False):
        for tp in timepoints_labels:
            for gr in groups_labels:
                # Only generate targets if there are eligible MAGs
                # Use the unpaired eligibility type for LMM
                mags = get_mags_by_eligibility(tp, gr, eligibility_type="lmm")
                if mags:  # Only proceed if there are eligible MAGs
                    group_str = ""  # no group marker for LMM
                    for taxon in tax_levels:
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "scores",
                                "processed",
                                "combined",
                                f"scores_lmm-{tp}-{gr}{group_str}-{taxon}.tsv",
                            )
                        )
                    
    return targets


def get_outlier_gene_targets():
    targets = []
    
    # Add significance test outlier targets if enabled
    if config["analysis"].get("use_significance_tests", True):
        for tp in timepoints_labels:
            for gr in groups_labels:
                for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                    group_str = ""  # no group marker for two-sample tests
                    mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                    for mag in mags:
                        prefix = f"{mag}_{test_type}{group_str}"
                        base_dir = os.path.join(
                            OUTDIR,
                            "outlier_genes",
                            f"{tp}-{gr}",
                        )
                        targets.append(
                            os.path.join(base_dir, f"{prefix}_outlier_genes.tsv")
                        )
                # Only include single sample targets if data_type is longitudinal
                if DATA_TYPE == "longitudinal":
                    sample_entries = get_single_sample_entries(tp, gr)
                    for mag, grp in sample_entries:
                        group_str = f"_{grp}"
                        prefix = f"{mag}_single_sample{group_str}"
                        base_dir = os.path.join(
                            OUTDIR,
                            "outlier_genes",
                            f"{tp}-{gr}",
                        )
                        targets.append(
                            os.path.join(base_dir, f"{prefix}_outlier_genes.tsv")
                        )
                        
    # Add LMM outlier targets if enabled
    if config["analysis"].get("use_lmm", False):
        for tp in timepoints_labels:
            for gr in groups_labels:
                group_str = ""  # no group marker for LMM
                # Use the unpaired eligibility type for LMM
                mags = get_mags_by_eligibility(tp, gr, eligibility_type="lmm")
                for mag in mags:
                    prefix = f"{mag}_lmm{group_str}"
                    base_dir = os.path.join(
                        OUTDIR,
                        "outlier_genes",
                        f"{tp}-{gr}",
                    )
                    targets.append(
                        os.path.join(base_dir, f"{prefix}_outlier_genes.tsv")
                    )
                    
    return targets


    
# Test the following functions before uncommenting. They are not all 100% up to date and might need modification. Be careful.
"""
def get_two_sample_targets(test_type):
    # Only generate targets if use_significance_tests is enabled
    if not config["analysis"].get("use_significance_tests", True):
        return []
        
    # Define subdirectory and file suffix based on the test type.
    if test_type == "two_sample_unpaired":
        subdir = "two_sample_unpaired"
        suffix = "_two_sample_unpaired.tsv.gz"
    elif test_type == "two_sample_paired":
        subdir = "two_sample_paired"
        suffix = "_two_sample_paired.tsv.gz"
    else:
        raise ValueError("test_type must be either 'unpaired' or 'paired'")
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag in get_mags_by_eligibility(tp, gr, eligibility_type=test_type):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "significance_tests",
                        f"{subdir}_{tp}-{gr}",
                        f"{mag}{suffix}",
                    )
                )
    return targets

def get_single_sample_targets():
    # Only generate targets if use_significance_tests is enabled
    if not config["analysis"].get("use_significance_tests", True):
        return []
        
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag, group in get_single_sample_entries(tp, gr):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "significance_tests",
                        f"single_sample_{tp}-{gr}",
                        f"{mag}_single_sample_{group}.tsv.gz",
                    )
                )
    return targets

def get_lmm_targets():
    # Only generate targets if use_lmm is enabled
    if not config["analysis"].get("use_lmm", False):
        return []
        
    targets = []
    for tp in timepoints_labels:
        for gr in groups_labels:
            for mag in get_mags_by_eligibility(tp, gr):
                targets.append(
                    os.path.join(
                        OUTDIR,
                        "significance_tests",
                        f"lmm_{tp}-{gr}",
                        f"{mag}_lmm.tsv.gz",
                    )
                )
    return targets

def get_significance_scores_targets():
    targets = []
    if config["analysis"].get("use_significance_tests", True):
        # For two-sample tests (unpaired and paired)
        for tp in timepoints_labels:
            for gr in groups_labels:
                for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                    for mag in get_mags_by_eligibility(tp, gr, eligibility_type=test_type):
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "scores",
                                "intermediate",
                                f"MAG_scores_{tp}-{gr}",
                                f"{mag}_score_{test_type}.tsv",
                            )
                        )
        # For single-sample tests - only include if data_type is longitudinal
        if DATA_TYPE == "longitudinal":
            for tp in timepoints_labels:
                for gr in groups_labels:
                    sample_entries = get_single_sample_entries(tp, gr)
                    for mag, group in sample_entries:
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "scores",
                                "intermediate",
                                f"MAG_scores_{tp}-{gr}",
                                f"{mag}_score_single_sample_{group}.tsv",
                            )
                        )
    # Add LMM scores
    if config["analysis"].get("use_lmm", False):
        for tp in timepoints_labels:
            for gr in groups_labels:
                group_str = ""  # no group marker for LMM
                # Use the unpaired eligibility type for LMM
                mags = get_mags_by_eligibility(tp, gr, eligibility_type="lmm")
                for mag in mags:
                    prefix = f"{mag}_lmm{group_str}"
                    base_dir = os.path.join(
                        OUTDIR,
                        "outlier_genes",
                        f"{tp}-{gr}",
                    )
                    targets.append(
                        os.path.join(base_dir, f"{prefix}_outlier_genes.tsv")
                    )
    return targets
    
def get_combined_scores_targets():
    targets = []
    # Two-sample tests: group_str is empty.
    if config["analysis"].get("use_significance_tests", True):
        for tp in timepoints_labels:
            for gr in groups_labels:
                for test_type in ["two_sample_unpaired", "two_sample_paired"]:
                    # Only generate targets if there are eligible MAGs for this test type
                    mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                    if mags:  # Only proceed if there are eligible MAGs
                        group_str = ""  # no group marker for two-sample tests
                        targets.append(
                            os.path.join(
                                OUTDIR,
                                "scores",
                                "processed",
                                "combined",
                                f"scores_{test_type}-{tp}-{gr}{group_str}-MAGs.tsv",
                            )
                        )
        # Single-sample tests: group_str is "_" plus the sample group.
        # Only include if data_type is longitudinal
        if DATA_TYPE == "longitudinal":
            for tp in timepoints_labels:
                for gr in groups_labels:
                    # get_single_sample_entries returns (mag, group) pairs for a given timepoint and group.
                    sample_entries = get_single_sample_entries(tp, gr)
                    if sample_entries:  # Only proceed if there are eligible MAGs
                        unique_groups = sorted(set([grp for mag, grp in sample_entries]))
                        for grp in unique_groups:
                            group_str = f"_{grp}"
                            targets.append(
                                os.path.join(
                                    OUTDIR,
                                    "scores",
                                    "processed",
                                    "combined",
                                    f"scores_single_sample-{tp}-{gr}{group_str}-MAGs.tsv",
                                )
                            )
    # Add LMM targets if enabled
    if config["analysis"].get("use_lmm", False):
        for tp in timepoints_labels:
            for gr in groups_labels:
                # Only generate targets if there are eligible MAGs
                # Use the unpaired eligibility type for LMM
                mags = get_mags_by_eligibility(tp, gr, eligibility_type="lmm")
                if mags:  # Only proceed if there are eligible MAGs
                    group_str = ""  # no group marker for LMM
                    targets.append(
                        os.path.join(
                            OUTDIR,
                            "scores",
                            "processed",
                            "combined",
                            f"scores_lmm-{tp}-{gr}{group_str}-MAGs.tsv",
                        )
                    )
    return targets

def get_gene_scores_targets():
    targets = []
    if config["analysis"].get("use_significance_tests", True):
        for tp in timepoints_labels:
            for gr in groups_labels:
                for test_type in ["two_sample_unpaired", "two_sample_paired"]:

                    # Only generate targets if there are eligible MAGs
                    mags = get_mags_by_eligibility(tp, gr, eligibility_type=test_type)
                    if mags:  # Only proceed if there are eligible MAGs
                        group_str = ""  # no group marker for two-sample tests
                        for mag in mags:
                            prefix = f"{mag}_{test_type}{group_str}"
                            base_dir = os.path.join(
                                OUTDIR,
                                "scores",
                                "processed",
                                f"gene_scores_{tp}-{gr}",
                            )
                            targets.extend(
                                [
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_combined.tsv"
                                    ),
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_individual.tsv"
                                    ),
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_overlapping.tsv"
                                    ),
                                ]
                            )
                # Only include single sample targets if data_type is longitudinal
                if DATA_TYPE == "longitudinal":
                    sample_entries = get_single_sample_entries(tp, gr)
                    if sample_entries:  # Only proceed if there are eligible MAGs
                        for mag, grp in sample_entries:
                            group_str = f"_{grp}"
                            prefix = f"{mag}_single_sample{group_str}"
                            base_dir = os.path.join(
                                OUTDIR,
                                "scores",
                                "processed",
                                f"gene_scores_{tp}-{gr}",
                            )
                            targets.extend(
                                [
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_combined.tsv"
                                    ),
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_individual.tsv"
                                    ),
                                    os.path.join(
                                        base_dir, f"{prefix}_gene_scores_overlapping.tsv"
                                    ),
                                ]
                            )
    return targets
"""
