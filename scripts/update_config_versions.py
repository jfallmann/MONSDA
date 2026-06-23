#!/usr/bin/env python
"""Synchronize the "VERSION" field in MONSDA config files with the
versioneer-derived package version (single source of truth).

This mirrors what versioneer does for the Python package: instead of hard
coding a version, the value is derived from the current git tag/state so the
config files always fit the release they ship with.

MONSDA config files may contain ``#`` comments (i.e. they are not strict
JSON), therefore replacement is line/regex based and preserves comments,
indentation and trailing commas. Only values that are empty, ``FIXME`` or
version-like are replaced, so documentation strings such as
``"VERSION": "Version of MONSDA to use, ..."`` are left untouched.
"""

import argparse
import subprocess
import sys
from pathlib import Path
import re

VERSION_LINE = re.compile(
    r'^(?P<pre>\s*"VERSION"\s*:\s*")(?P<val>[^"]*)(?P<post>".*)$', re.MULTILINE
)
VERSION_LIKE = re.compile(r"^v?\d+(\.\d+)*([._\-+].*)?$")


def get_version():
    """Return the versioneer-derived version (same value used at runtime)."""
    try:
        from MONSDA import __version__

        return __version__
    except Exception:
        root = Path(__file__).resolve().parent.parent
        out = subprocess.check_output(
            [sys.executable, "setup.py", "--version"], cwd=root
        )
        return out.decode().strip().splitlines()[-1].strip()


def is_replaceable(val):
    v = val.strip()
    return v == "" or v.upper() == "FIXME" or bool(VERSION_LIKE.match(v))


def update_text(text, version):
    changed = {"n": 0}

    def repl(m):
        if is_replaceable(m.group("val")) and m.group("val") != version:
            changed["n"] += 1
            return m.group("pre") + version + m.group("post")
        return m.group(0)

    return VERSION_LINE.sub(repl, text), changed["n"]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", nargs="+", help="Config files to update")
    parser.add_argument(
        "--version",
        default=None,
        help="Version to write (default: versioneer-derived package version)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Do not modify files; exit non-zero if any file is out of sync",
    )
    args = parser.parse_args(argv)

    version = args.version or get_version()
    rc = 0
    for f in args.files:
        path = Path(f)
        if not path.is_file():
            print(f"skip (not found): {f}", file=sys.stderr)
            continue
        text = path.read_text()
        new, n = update_text(text, version)
        if n == 0:
            continue
        if args.check:
            print(f"out of sync: {f} (expected VERSION={version})", file=sys.stderr)
            rc = 1
        else:
            path.write_text(new)
            print(f"updated {f}: VERSION -> {version} ({n} field(s))")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
