import os
from pathlib import Path

import nox


def get_python_version():
    """Get the python version from .python-version file unless it is github actions."""
    if "GITHUB_ACTIONS" in os.environ:
        return Path(os.environ["Python_ROOT_DIR"]).parent.name
    return Path(".python-version").read_text().strip()


def get_source_dirs():
    """Every folder holding pipeline code: the numbered steps, the alternative pipeline, species_specific/, scripts/ and src/.

    There is no single package here, so this replaces the PACKAGE_NAME target. results/ only holds
    snapshot copies of pipeline scripts.
    """
    dirs = sorted(p for p in Path(".").glob("[0-9][0-9]_*") if p.is_dir())
    dirs += [Path("Alternative_Get_PanLoci_pipeline/scripts"), Path("species_specific"), Path("scripts"), Path("src")]
    return [str(d) for d in dirs if d.is_dir()]


def get_script_dirs():
    """mypy targets: the src/ package as one, plus every other folder that directly contains .py files.

    Loose scripts outside the package are checked one folder at a time, as some basenames repeat.
    """
    loose = {str(p.parent) for d in SOURCE_DIRS if d != "src" for p in Path(d).rglob("*.py")}
    return ["src", *sorted(loose)]


def get_toml_files():
    return [str(p) for p in Path(".").rglob("*.toml") if ".venv" not in p.parts and ".nox" not in p.parts]


PYTHON_VERSION = get_python_version()
SOURCE_DIRS = get_source_dirs()

nox.options.reuse_existing_virtualenvs = True
nox.options.default_venv_backend = "uv"


def uv(session, *args):
    """Run uv commands inside a nox session."""
    session.run("uv", *args, "--active", external=True)


@nox.session(python=PYTHON_VERSION)
def install(session):
    """
    Equivalent of:
      make install
    """
    uv(session, "sync", "--group", "dev")


@nox.session(python=PYTHON_VERSION)
def lock(session):
    """
    Equivalent of:
      make lock
    """
    # uv lock doesn't accept --active, so it's called directly rather than via uv().
    session.run("uv", "lock", "--upgrade", external=True)


@nox.session(python=PYTHON_VERSION)
def lock_check(session):
    """
    Equivalent of:
      make lock-check
    """
    session.run("uv", "lock", "--check", external=True)


#
# Formatting / chores
#


@nox.session(python=PYTHON_VERSION)
def ruff_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "check", ".", "--fix")


@nox.session(python=PYTHON_VERSION)
def black_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "format", ".")


@nox.session(python=PYTHON_VERSION)
def tomlsort_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("tombi", "format", *get_toml_files(), external=True)


@nox.session(python=PYTHON_VERSION)
def isort_fixes(session):
    uv(session, "sync", "--group", "dev")
    session.run("isort", *SOURCE_DIRS, "tests", "noxfile.py")


@nox.session(python=PYTHON_VERSION)
def chores(session):
    """
    Equivalent of:
      make chores

    This runs all formatting and linting fixes but also mypy checks.
    """
    session.notify("isort_fixes")
    session.notify("ruff_fixes")
    session.notify("black_fixes")
    session.notify("tomlsort_fixes")
    session.notify("mypy_check")


#
# Checks / linting
#


@nox.session(python=PYTHON_VERSION)
def ruff_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "check")


@nox.session(python=PYTHON_VERSION)
def black_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("ruff", "format", ".", "--check")


@nox.session(python=PYTHON_VERSION)
def mypy_check(session):
    uv(session, "sync", "--group", "dev")
    failed = []
    for script_dir in get_script_dirs():
        try:
            session.run("mypy", script_dir)
        except nox.command.CommandFailed:
            failed.append(script_dir)
    if failed:
        session.error("mypy failed in: " + ", ".join(failed))


@nox.session(python=PYTHON_VERSION)
def tomlsort_check(session):
    uv(session, "sync", "--group", "dev")
    tomls = get_toml_files()
    session.run("tombi", "lint", *tomls, external=True)
    session.run("tombi", "format", *tomls, "--check", external=True)


@nox.session(python=PYTHON_VERSION)
def isort_check(session):
    uv(session, "sync", "--group", "dev")
    session.run("isort", *SOURCE_DIRS, "tests", "noxfile.py", "--check-only")


#
# Testing
#


@nox.session(python=PYTHON_VERSION)
def tests(session):
    """
    Equivalent of:
      make tests
    """
    session.notify("ruff_check")
    session.notify("black_check")
    session.notify("mypy_check")
    session.notify("tomlsort_check")
    session.notify("isort_check")
    session.notify("pytest")


@nox.session(python=PYTHON_VERSION)
def precommit(session):
    """
    Equivalent of:
      make precommit
    """
    session.notify("ruff_check")
    session.notify("black_check")
    session.notify("mypy_check")
    session.notify("tomlsort_check")
    session.notify("isort_check")


@nox.session(python=PYTHON_VERSION)
def pytest(session):
    uv(session, "sync", "--group", "dev")
    session.run(
        "pytest",
        "--cov=.",
        "--cov-report=term-missing",
        "tests",
    )


@nox.session(python=PYTHON_VERSION)
def pytest_loud(session):
    uv(session, "sync", "--group", "dev")
    session.run(
        "pytest",
        "--log-cli-level=DEBUG",
        "-o",
        "log_cli=true",
        "--cov=.",
        "--cov-report=term-missing",
        "tests",
    )
