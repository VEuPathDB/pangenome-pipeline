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
