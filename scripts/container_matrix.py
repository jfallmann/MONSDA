#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from MONSDA.Containers import normalize_container_version  # noqa: E402


def environment_names(root):
    names = sorted(path.stem for path in (root / "envs").glob("*.yaml"))
    if not names:
        raise RuntimeError("No environment definitions found")
    return names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--version", action="store_true")
    args = parser.parse_args()
    if args.version:
        print(f"version={normalize_container_version()}")
    else:
        print(json.dumps({"environment": environment_names(args.root)}))


if __name__ == "__main__":
    main()
