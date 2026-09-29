"""step04_interproscan/update_interpro_in_gff.py - replacing domain Dbxrefs."""

import pytest

SCRIPT = "src/pangenome_pipeline/step04_interproscan/update_interpro_in_gff.py"


@pytest.fixture(scope="module")
def ipr(load_script):
    return load_script(SCRIPT)


def test_strip_old_domain_xrefs_keeps_other_databases(ipr):
    assert ipr.strip_old_domain_xrefs("GeneID:1,InterPro:IPR000001,pfam:PF00001,UniProt:Q1") == [
        "GeneID:1",
        "UniProt:Q1",
    ]


def test_new_domain_tokens(ipr):
    assert ipr.new_domain_tokens("IPR1;IPR2", "PF1") == ["InterPro:IPR1", "InterPro:IPR2", "PFAM:PF1"]
    assert ipr.new_domain_tokens("", "") == []


def test_load_lookup_skips_header_and_short_rows(ipr, write_file):
    path = write_file("lookup.tsv", "protein\tinterpro\tpfam\np1\tIPR1\tPF1\nbad_row\n")
    assert ipr.load_lookup(path) == {"p1": ("IPR1", "PF1")}


def test_manifest_job_reads_and_writes_step04_out(run_job, tmp_path):
    """Submitted from a project folder: proteins and manifest both live in output/step04_interproscan/."""
    project = tmp_path / "project"
    proteins = project / "output" / "step04_interproscan" / "proteins"
    proteins.mkdir(parents=True)
    (project / "genome_order.txt").write_text("G1\nG2\n")
    (proteins / "G1.faa").write_text("".join(f">p{i}\nM\n" for i in range(3)))
    (proteins / "G2.faa").write_text(">q1\nM\n")

    result = run_job("workflow/step04_interproscan/build_interpro_manifest.slurm.sh", project)
    assert result.returncode == 0, result.stderr
    manifest = (project / "output" / "step04_interproscan" / "interpro_manifest.tsv").read_text().splitlines()
    assert manifest[1:] == ["1\tG1\t1\t1\t1000", "2\tG2\t1\t1\t1000"]
    assert "Total chunks: 2" in result.stdout
