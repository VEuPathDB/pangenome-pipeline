"""Helpers shared by the pipeline's step scripts (standard library only).

The step scripts aren't installed as a package, so each one that uses these adds the repo's
lib/ folder to sys.path first, relative to its own location:

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
    from pangenome_common.intervals import merged_len
"""
