"""Reading the small list/config files the pipeline passes between steps."""

import os


def read_strains(value):
    """Genome names, in order: from a file with one per line ('#' lines and blank lines skipped),
    or from a comma-separated string if `value` isn't a file. The first one is the reference."""
    if os.path.isfile(value):
        with open(value) as fh:
            return [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    return [s.strip() for s in value.split(",") if s.strip()]
