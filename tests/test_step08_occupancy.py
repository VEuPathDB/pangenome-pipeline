"""08_Pangenome_Visualisation - pangene occupancy categories and matrix reading."""

import pytest

SCRIPTS = "08_Pangenome_Visualisation/scripts/"


@pytest.fixture(scope="module")
def occ(load_script):
    return load_script(SCRIPTS + "plot_occupancy.py")


@pytest.fixture(scope="module")
def pancore(load_script):
    return load_script(SCRIPTS + "plot_pancore.py")


@pytest.mark.parametrize(
    ("present_in", "expected"),
    [(10, "Core"), (9, "Soft-core"), (8, "Shell"), (3, "Shell"), (2, "Cloud"), (1, "Cloud")],
)
def test_categorize_defaults(occ, present_in, expected):
    assert occ.categorize(present_in, 10) == expected


def test_categorize_rejects_inverted_thresholds(occ):
    with pytest.raises(ValueError):
        occ.categorize(5, 10, core_missing=1, softcore_missing=1)


def test_read_presence(pancore, write_file):
    text = (
        "cluster\tA_for_get_pangenes\tB_for_get_pangenes\tC\n"
        "c1\tgeneA1\tgeneB1\t-\n"
        "c2\t-\tNA\t\n"  # absent everywhere -> dropped
        "c3\t0\tgeneB3\tgeneC3\n"
    )
    genomes, pres = pancore.read_presence(write_file("matrix.tab", text))
    assert genomes == ["A", "B", "C"]
    assert pres.tolist() == [[True, True, False], [False, True, True]]
