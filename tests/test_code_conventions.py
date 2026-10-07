"""Code conventions ruff has no rule for.

Constants at the top: module-level UPPER_CASE names (thresholds, regexes, lookup tables, colours)
are defined before the first function or class, so a script's settings are all in one place.
A constant computed by calling one of the module's own functions (e.g. noxfile.py's
PYTHON_VERSION = get_python_version()) can't move above that function and is allowed - as is a
constant built from such a constant.
"""

import ast
import re
import subprocess

from conftest import REPO_ROOT

CONSTANT_NAME = re.compile(r"^_?[A-Z][A-Z0-9_]*$")


def _python_files():
    files = subprocess.run(
        ["git", "ls-files", "*.py"], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    return [f for f in files if not f.startswith("results/")]


def _assigned_names(node):
    if isinstance(node, ast.Assign):
        targets = node.targets
    elif isinstance(node, ast.AnnAssign):
        targets = [node.target]
    else:
        return []
    names = []
    for target in targets:
        elts = target.elts if isinstance(target, ast.Tuple) else [target]
        names += [t.id for t in elts if isinstance(t, ast.Name)]
    return names


def misplaced_constants(source):
    """(line, name) for module-level constants defined after the first def/class."""
    tree = ast.parse(source)
    defs = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))]
    if not defs:
        return []
    first_def = defs[0].lineno
    must_follow = {n.name for n in defs}  # own functions, then constants computed from them
    found = []
    for node in tree.body:
        if node.lineno < first_def:
            continue
        value = getattr(node, "value", None)
        needs_a_def = value is not None and any(
            isinstance(n, ast.Name) and n.id in must_follow for n in ast.walk(value)
        )
        for name in _assigned_names(node):
            if needs_a_def:
                must_follow.add(name)
            elif CONSTANT_NAME.match(name):
                found.append((node.lineno, name))
    return found


def test_constants_are_defined_before_functions():
    hits = [
        f"{f}:{line}: {name}"
        for f in _python_files()
        for line, name in misplaced_constants((REPO_ROOT / f).read_text())
    ]
    assert not hits, "Move these constants above the first def/class:\n  " + "\n  ".join(hits)


def test_misplaced_constants_detection():
    source = (
        "import re\n"
        "LIMIT = 5\n"
        "def f():\n    return 1\n"
        "COLOURS = {'a': 1}\n"
        "PATTERN: re.Pattern = re.compile('x')\n"
        "COMPUTED = f()\n"  # needs f - allowed
        "DERIVED = COMPUTED + 1\n"  # built from COMPUTED - allowed
        "lower_case = 2\n"  # not a constant
        "A, B = 1, 2\n"
    )
    assert misplaced_constants(source) == [(5, "COLOURS"), (6, "PATTERN"), (10, "A"), (10, "B")]
