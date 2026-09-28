"""step08_pangenome_visualisation - pangene occupancy categories, matrix reading and plot_all.sh."""

import os
import subprocess
import sys

import pytest
from conftest import REPO_ROOT

SCRIPTS = "src/pangenome_pipeline/step08_pangenome_visualisation/"


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


def test_plot_all_runs_every_plot(tmp_path):
    """plot_all.sh, run directly, finds the package itself and writes every figure."""
    rows = ["cluster\tA_for_get_pangenes\tB_for_get_pangenes\tC_for_get_pangenes\tD_for_get_pangenes"]
    for i in range(40):  # a spread of occupancies: core, shell and cloud clusters
        present = [(i + j) % 4 != 0 or i % 3 == 0 for j in range(4)]
        rows.append(f"c{i}\t" + "\t".join(f"g{i}_{j}" if p else "-" for j, p in enumerate(present)))
    matrix = tmp_path / "pangene_matrix.tr.tab"
    matrix.write_text("\n".join(rows) + "\n")
    venv_bin = os.path.dirname(sys.executable)  # stands in for the environment the job activates (matplotlib, scipy)
    result = subprocess.run(
        [
            "bash",
            str(REPO_ROOT / "08_Pangenome_Visualisation" / "scripts" / "plot_all.sh"),
            "-i",
            str(matrix),
            "-o",
            str(tmp_path / "out"),
            "-l",
            "T",
        ],
        cwd=tmp_path,
        env={**os.environ, "MPLBACKEND": "Agg", "PATH": f"{venv_bin}:{os.environ['PATH']}"},
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    for name in ("T_pancore.png", "T_pancore_separate.png", "T_occupancy.png", "T_singletons.png"):
        assert (tmp_path / "out" / name).is_file(), name
