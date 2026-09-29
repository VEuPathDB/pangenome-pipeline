"""docs/migrating_projects.md's move script, run on a project in the old layout."""

import re
import subprocess

from conftest import REPO_ROOT

OLD_LAYOUT = [
    "species_check/species_check.tsv",
    "01_Genome_Selection_QC/contig_size_profile/contig_sizes.tsv",
    "01_Genome_Selection_QC/drop_contigs.txt",
    "01_Genome_Selection_QC/scripts/.keep",
    "RM_AF293/AF293.fna.out",
    "curated_TE_library.fa",
    "02_TE_annotation_Clean_Up/final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa",
    "03_denovo_Genome_Annotation/funannotate/GENOME.sorted.fna",
    "03_denovo_Genome_Annotation/README.md",
    "interpro/interpro_manifest.tsv",
    "gff_with_interpro/AF293.gff3",
    "liftoff_output/polished/a_onto_b_liftoff.gff3_polished",
    "filtered_annotation/AF293_final_annotation_filtered.gff",
    "for_get_pangenes_7b_pangenes/x/pangene_matrix.tr.tab",
    "pangenome_plots/Afum7b_pancore.png",
    "genomes/AF293.fna",
]
NEW_LAYOUT = [
    "output/step01_genome_selection_qc/species_check/species_check.tsv",
    "output/step01_genome_selection_qc/contig_size_profile/contig_sizes.tsv",
    "output/step01_genome_selection_qc/drop_contigs.txt",
    "01_Genome_Selection_QC/scripts/.keep",  # not a result: left alone
    "output/step02_te_annotation_clean_up/repeatmasker/RM_AF293/AF293.fna.out",
    "output/step02_te_annotation_clean_up/curated_TE_library.fa",
    "output/step02_te_annotation_clean_up/te_library/final_TE_candidates_after_MCHelper/final_potential_new_TEs.fa",
    "output/step03_denovo_genome_annotation/funannotate/GENOME.sorted.fna",
    "03_denovo_Genome_Annotation/README.md",
    "output/step04_interproscan/interpro_manifest.tsv",
    "output/step04_interproscan/gff_with_interpro/AF293.gff3",
    "output/step05_cross_genome_annotation_integration/liftoff_output/polished/a_onto_b_liftoff.gff3_polished",
    "output/step05_cross_genome_annotation_integration/filtered_annotation/AF293_final_annotation_filtered.gff",
    "output/step06_get_pangenes/for_get_pangenes_7b_pangenes/x/pangene_matrix.tr.tab",
    "output/step08_pangenome_visualisation/pangenome_plots/Afum7b_pancore.png",
    "genomes/AF293.fna",  # inputs stay
]


def _script():
    text = (REPO_ROOT / "docs" / "migrating_projects.md").read_text()
    section = text[text.index("## 2. Move the results") :]
    return re.search(r"```bash\n(.*?)```", section, re.DOTALL).group(1)


def test_move_script_moves_every_result_and_is_safe_to_rerun(tmp_path):
    for rel in OLD_LAYOUT:
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text(rel)
    # something already at a new place is never overwritten
    kept = tmp_path / "output/step05_cross_genome_annotation_integration/filtered_annotation"
    kept.mkdir(parents=True)
    (kept / "newer.gff").write_text("newer")

    run = subprocess.run(
        ["bash", "-c", "set -u\n" + _script()], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert run.returncode == 0, run.stderr
    for old, new in zip(OLD_LAYOUT, NEW_LAYOUT, strict=True):
        if "filtered_annotation" in old:
            assert (tmp_path / old).exists() and "SKIP  filtered_annotation" in run.stdout
            continue
        assert (tmp_path / new).read_text() == old, new
    assert (kept / "newer.gff").read_text() == "newer"
    assert (tmp_path / "logs").is_dir()

    again = subprocess.run(
        ["bash", "-c", "set -u\n" + _script()], cwd=tmp_path, capture_output=True, text=True, check=False
    )
    assert again.returncode == 0, again.stderr
    assert "MOVED" not in again.stdout
