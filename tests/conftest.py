"""Shared test helpers.

The pipeline scripts aren't a package - they live in step folders, some with hyphens in the name,
and import their siblings via sys.path. load_script() imports one by file path the same way.
"""

import importlib.util
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("MPLBACKEND", "Agg")

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load(relpath):
    path = REPO_ROOT / relpath
    script_dir = str(path.parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)  # sibling imports, e.g. `from plot_pancore import read_presence`
    name = "pipeline_" + path.stem
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def load_script():
    return _load


@pytest.fixture
def write_file(tmp_path):
    """Write text to tmp_path/<name> and return the path."""

    def _write(name, text):
        path = tmp_path / name
        path.write_text(text)
        return path

    return _write


def gff_line(seqid, ftype, start, end, strand, attrs, source="test"):
    return "\t".join([seqid, source, ftype, str(start), str(end), ".", strand, ".", attrs]) + "\n"
