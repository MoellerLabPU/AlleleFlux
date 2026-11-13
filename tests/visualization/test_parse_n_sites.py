"""``parse_n_sites``: the "how many top sites" option of alleleflux-plot-trajectories.

Valid: the word ``all`` in any letter case (every site), or a positive whole
number (the top N).  Anything else raises, so a typo in the config fails the
plotting job loudly instead of exiting 0 with no output (which Snakemake then
reports as a missing output after a job that looked successful; found 2026-10-19
when the rule accepted ``ALL`` but the script only accepted ``all``).
"""

import pytest

from alleleflux.scripts.visualization.plot_allele_trajectory import parse_n_sites


@pytest.mark.parametrize("word", ["all", "ALL", "All", " all "])
def test_the_word_all_in_any_case_means_every_site(word):
    assert parse_n_sites(word, "n_sites_line") == "all"


@pytest.mark.parametrize(
    "value, expected", [("2000", 2000), ("1", 1), (7, 7), (" 15 ", 15)]
)
def test_a_positive_whole_number_is_the_top_n(value, expected):
    assert parse_n_sites(value, "n_sites_line") == expected


@pytest.mark.parametrize("bad", ["0", "-3", "2.5", "alll", "", "none", "1e3"])
def test_anything_else_raises_naming_the_option(bad):
    with pytest.raises(ValueError, match="n_sites_dist"):
        parse_n_sites(bad, "n_sites_dist")
