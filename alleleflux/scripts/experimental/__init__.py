# Experimental subpackage for AlleleFlux
#
# Modules here are NOT used or benchmarked in the AlleleFlux publication and are
# disabled by default in the shipped config template. They are kept because they
# may be useful later, but their results have not been validated -- treat any
# output as provisional. A module graduates out of this folder (back into
# analysis/, preprocessing/, ...) once it has been validated and tested.
#
# Current contents:
#   regional_contrast.py          -- gene/window-level contrast of allele-frequency
#                                    change between groups (alleleflux-regional-contrast).
#                                    Suspected to not behave as intended; needs review
#                                    before any use.
#   regional_contrast_summary.py  -- summarises regional_contrast output
#                                    (alleleflux-regional-contrast-summary).
#   outliers_genes.py             -- outlier gene detection via binomial and Poisson
#                                    tests of gene vs MAG scores (alleleflux-outliers).
