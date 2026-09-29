"""alternative_get_panloci/panloci_common.py - helpers shared by step1-3."""

import os
import subprocess

import pytest
from conftest import REPO_ROOT

SCRIPT = "src/pangenome_pipeline/alternative_get_panloci/panloci_common.py"

SAMPLES = (
    "strain_id\tgenome_fasta\tannotation_gff\tanalysis_order\nB \tb.fa\tNA\t2\nA\ta.fa\tgffs/\t1\nC\tc.fa\tNA\t10\n"
)


@pytest.fixture(scope="module")
def pc(load_script):
    return load_script(SCRIPT)


def test_parse_samples_sorted_or_file_order(pc, write_file):
    path = write_file("samples.tsv", SAMPLES)
    assert [s["strain_id"] for s in pc.parse_samples(path)] == ["A", "B", "C"]
    assert [s["strain_id"] for s in pc.parse_samples(path, sort_by_order=False)] == ["B", "A", "C"]
    assert pc.parse_samples(path)[0] == {
        "strain_id": "A",
        "genome_fasta": "a.fa",
        "annotation_gff": "gffs/",
        "analysis_order": 1,
    }


def test_parse_samples_errors(pc, write_file, tmp_path):
    with pytest.raises(FileNotFoundError):
        pc.parse_samples(str(tmp_path / "missing.tsv"))
    with pytest.raises(ValueError, match="analysis_order"):
        pc.parse_samples(write_file("bad.tsv", "strain_id\tgenome_fasta\tannotation_gff\nA\ta.fa\tNA\n"))


def test_find_gff_file(pc, tmp_path):
    for name in ("3D7_annot.gff3", "DD2.gff", "other.gff"):
        (tmp_path / name).write_text("")
    d = str(tmp_path)
    assert pc.find_gff_file(d, "REF3D7").endswith("3D7_annot.gff3")  # REF prefix dropped
    assert pc.find_gff_file(d, "dd2").endswith("DD2.gff")  # upper-case variant
    assert pc.find_gff_file(d, "HB3") is None  # several GFFs, none matching
    assert pc.find_gff_file(d, "REF") is None  # "REF" minus "REF" is "", which must not match everything
    assert pc.find_gff_file("NA", "A") is None
    assert pc.find_gff_file(str(tmp_path / "DD2.gff"), "whatever").endswith("DD2.gff")


def test_find_gff_file_single_candidate_fallback(pc, tmp_path):
    (tmp_path / "only.gff").write_text("")
    assert pc.find_gff_file(str(tmp_path), "XYZ").endswith("only.gff")


def test_load_fasta_and_missing(pc, write_file, tmp_path, capsys):
    path = write_file("g.fa", ">chr1 desc\nacgT\nNNa\n\n>chr2\nGG\n")
    assert pc.load_fasta(path) == {"chr1": "ACGTNNA", "chr2": "GG"}
    assert pc.load_fasta(str(tmp_path / "none.fa"), missing_ok=True) is None
    assert "FASTA file not found" in capsys.readouterr().err
    with pytest.raises(FileNotFoundError):
        pc.load_fasta(str(tmp_path / "none.fa"))


def test_extract_sequence(pc):
    fasta = {"chr1": "ACGTNNA"}
    assert pc.extract_sequence(fasta, "chr1", 0, 3) == "ACG"
    assert pc.extract_sequence(fasta, "chr1", 2, 99, "-") == "TNNAC"  # clipped, reverse-complemented
    assert pc.extract_sequence(fasta, "chr1", 5, 5) is None
    assert pc.extract_sequence(fasta, "chrX", 0, 1) is None
    assert pc.extract_sequence(None, "chr1", 0, 1) is None


def test_intervals_and_attributes(pc):
    assert pc.merge_overlapping_intervals([(20, 30), (0, 10), (10, 15), (5, 8)]) == [(0, 15), (20, 30)]
    assert pc.parse_gff_attributes(" ID=a;;Parent=b,c ") == {"ID": "a", "Parent": "b,c"}
    assert pc.split_parent_ids("b,,c") == ["b", "c"]
    assert pc.is_pseudogene_feature("gene", {"biotype": "Processed_Pseudogene"})
    assert not pc.is_pseudogene_feature("gene", {"gene_biotype": "protein_coding"})


def test_cluster_mapped_blocks_by_gap_and_strand(pc):
    rows = [
        ["chr1", "100", "200", "g", "0", "+"],
        ["chr1", "250", "300", "g", "0", "+"],  # 50 bp gap
        ["chr1", "1000", "1100", "g", "0", "+"],  # 700 bp gap
        ["chr1", "260", "280", "g", "0", "-"],
    ]
    clusters = sorted(pc.cluster_mapped_blocks(rows, max_gap=100))
    assert clusters == [
        (("chr1", "+"), [(100, 200), (250, 300)]),
        (("chr1", "+"), [(1000, 1100)]),
        (("chr1", "-"), [(260, 280)]),
    ]


def test_merge_hal_liftover_primary_and_duplicates(pc, write_file, tmp_path):
    seed = write_file("seed.bed", "chr1\t100\t400\tgA\t0\t+\nchr1\t500\t600\tgA\t0\t+\n")
    raw = write_file(
        "raw.bed",
        "chr1\t1000\t1200\tgA\t0\t+\nchr1\t1150\t1300\tgA\t0\t+\nchr5\t10\t60\tgA\t0\t+\nshort\n",
    )
    out, qc = tmp_path / "merged.bed", tmp_path / "multihit.tsv"
    assert pc.merge_hal_liftover(raw, seed, str(out), str(qc), 5000, id_column="panlocus_id") == (1, 1)
    # primary = the chr1 cluster: 300 merged bp of a 400 bp seed
    assert out.read_text() == "chr1\t1000\t1300\tgA\t0\t+\t400\t300\t75.0\n"
    assert qc.read_text().splitlines() == [
        "panlocus_id\tdup_chrom\tdup_start\tdup_end\tdup_strand\torig_size\tdup_mapped_len\tdup_percentage\thit_rank",
        "gA\tchr5\t10\t60\t+\t400\t50\t12.5\t2",
    ]
    pc.merge_hal_liftover(raw, seed, str(out), str(qc), 5000)  # appends, no second header
    assert len(qc.read_text().splitlines()) == 3


def test_runner_finds_the_pipeline_and_defaults_to_alt_panloci_out(tmp_path):
    """run_get_panloci.sh finds the checkout from its own folder; with no RUN_DIR set it works in
    output/alternative_get_panloci/ under the folder it is run from (and stops there: no samples.tsv)."""
    runner = REPO_ROOT / "workflow" / "alternative_get_panloci" / "no_slurm" / "run_get_panloci.sh"
    env = {k: v for k, v in os.environ.items() if k not in ("PIPELINE_DIR", "RUN_DIR", "SLURM_SUBMIT_DIR")}
    env["PIPELINE_CONFIG"] = str(tmp_path / "no_config.sh")
    found = subprocess.run(["bash", str(runner)], cwd=tmp_path, env=env, capture_output=True, text=True, check=False)
    assert found.returncode == 1
    assert f"Samples table not found: {tmp_path}/output/alternative_get_panloci/samples.tsv" in found.stderr

    wrong = subprocess.run(
        ["bash", str(runner)],
        cwd=tmp_path,
        env={**env, "PIPELINE_DIR": str(tmp_path)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert wrong.returncode == 1
    assert "pipeline not found - export PIPELINE_DIR" in wrong.stderr
