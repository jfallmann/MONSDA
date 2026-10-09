#!/usr/bin/env python3
"""Dump shell/script blocks from workflows/*.smk and workflows/*.nf.

Reusable version of the throwaway extractor used to build
smk_shell_commands.txt / nf_script_commands.txt during the manual CLI-flag
audit. Regenerate these dumps with --refresh whenever a workflow's shell
commands change, then use them as a reference while updating
tool_manifest.yaml by hand (the dumps themselves are not consumed by
audit_tool_flags.py).

Usage:
    python tests/tool_audit/extract_workflow_commands.py --refresh
    python tests/tool_audit/extract_workflow_commands.py --check-manifest
"""

import argparse
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
WORKFLOWS_DIR = REPO_ROOT / "workflows"

SMK_SHELL_RE_DOUBLE = re.compile(r'shell:\s*"((?:[^"\\]|\\.)*)"')
SMK_SHELL_RE_TRIPLE = re.compile(r"shell:\s*'''(.*?)'''", re.S)
NF_SCRIPT_RE = re.compile(r'"""(.*?)"""', re.S)


def iter_shell_blocks(smk_path):
    text = smk_path.read_text()
    for m in SMK_SHELL_RE_DOUBLE.finditer(text):
        yield m.group(1)
    for m in SMK_SHELL_RE_TRIPLE.finditer(text):
        yield m.group(1)


def iter_script_blocks(nf_path):
    text = nf_path.read_text()
    for m in NF_SCRIPT_RE.finditer(text):
        yield m.group(1)


def dump_smk(out_path):
    with open(out_path, "w") as out:
        for f in sorted(WORKFLOWS_DIR.glob("*.smk")):
            for block in iter_shell_blocks(f):
                out.write(f"### {f.name}\n{block}\n\n")


def dump_nf(out_path):
    with open(out_path, "w") as out:
        for f in sorted(WORKFLOWS_DIR.glob("*.nf")):
            for block in iter_script_blocks(f):
                out.write(f"### {f.name}\n{block}\n\n")


def check_manifest_coverage(manifest_path, smk_dump, nf_dump):
    """Sanity check: every source_ref in the manifest should still be found
    as a filename header in the corresponding dump, catching the extractor
    silently failing to parse new syntax."""
    import yaml

    with open(manifest_path) as fh:
        manifest = yaml.safe_load(fh).get("tools", {})

    smk_text = smk_dump.read_text() if smk_dump.exists() else ""
    nf_text = nf_dump.read_text() if nf_dump.exists() else ""

    problems = []
    for tool, spec in manifest.items():
        for inv in spec.get("invocations", []):
            for ref in inv.get("source_refs", []):
                name = Path(ref).name
                if ref.endswith(".smk") and f"### {name}" not in smk_text:
                    problems.append(f"{tool}: {ref} not found in {smk_dump.name}")
                elif ref.endswith(".nf") and f"### {name}" not in nf_text:
                    problems.append(f"{tool}: {ref} not found in {nf_dump.name}")
    return problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--refresh", action="store_true",
                     help="regenerate smk_shell_commands.txt / nf_script_commands.txt")
    ap.add_argument("--check-manifest", action="store_true",
                     help="check tool_manifest.yaml source_refs are still covered by the dumps")
    args = ap.parse_args()

    smk_dump = HERE / "smk_shell_commands.txt"
    nf_dump = HERE / "nf_script_commands.txt"

    if args.refresh:
        dump_smk(smk_dump)
        dump_nf(nf_dump)
        print(f"wrote {smk_dump}")
        print(f"wrote {nf_dump}")

    if args.check_manifest:
        problems = check_manifest_coverage(HERE / "tool_manifest.yaml", smk_dump, nf_dump)
        if problems:
            print("Manifest coverage problems:")
            for p in problems:
                print(f"  - {p}")
            return 1
        print("Manifest source_refs all covered by current dumps.")

    if not args.refresh and not args.check_manifest:
        ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
