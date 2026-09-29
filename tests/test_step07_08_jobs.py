"""Steps 07 and 08 read get_pangenes' output from step 06's output folder by default."""

STEP06 = "output/step06_get_pangenes"


def test_plot_all_job_set_reads_step06_out(run_job, tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    result = run_job(
        "workflow/step08_pangenome_visualisation/plot_all_job.sh", project, "--set", "7b", "--label", "X", SKIP_ENV="1"
    )
    assert result.returncode == 1
    assert f"get_pangenes folder not found: {project}/{STEP06}/for_get_pangenes_7b_pangenes" in result.stderr


def test_explore_multigene_job_reads_step06_out(run_job, tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    result = run_job(
        "workflow/step07_post_get_pangenes_adjusting/explore_multigene_pangenes_job.sh", project, LABEL="X"
    )
    assert result.returncode == 1
    assert f"expected exactly one {project}/{STEP06}/for_get_pangenes_7b_pangenes/*2neigh_algMmap_/" in result.stderr
    assert (project / "output" / "step07_post_get_pangenes_adjusting" / "multigene_pangenes").is_dir()
