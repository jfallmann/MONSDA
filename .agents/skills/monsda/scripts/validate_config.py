#!/usr/bin/env python3
"""Cheap, offline sanity checks for a MONSDA config.json.

This does NOT replace MONSDA's own validation (``monsda --save`` + a
Snakemake/Nextflow dry run, see reference/cli.md) — it only catches obvious
mistakes (bad JSON, missing step blocks, a stale VERSION, sample/reference
files that don't exist) before you spend time on a real generate/run.

Usage:
    python validate_config.py CONFIG.json [--base DIR]

``--base`` is the run/project directory paths in the config are resolved
against (defaults to the config file's own directory).
"""
import argparse
import json
import os
import sys

KNOWN_WORKFLOW_STEPS = {
    "FETCH", "BASECALL", "QC", "TRIMMING", "MAPPING", "DEDUP", "COUNTING",
    "TRACKS", "PEAKS", "DE", "DEU", "DAS", "DTU", "CIRCS", "FUSIONS",
}


def fail(msg, errors):
    errors.append(msg)


def walk_settings_leaves(node, path=()):
    """Yield (path, leaf_dict) for every SETTINGS leaf that has a SAMPLES key."""
    if not isinstance(node, dict):
        return
    if "SAMPLES" in node:
        yield path, node
        return
    for key, val in node.items():
        if isinstance(val, dict):
            yield from walk_settings_leaves(val, path + (key,))


def check_installed_version(declared_version, errors, warnings):
    try:
        import MONSDA
        installed = getattr(MONSDA, "__version__", None)
    except Exception:
        warnings.append(
            "Could not import MONSDA to check VERSION against the installed "
            f"release (declared VERSION={declared_version!r}); verify manually "
            "with: python -c \"import MONSDA; print(MONSDA.__version__)\""
        )
        return
    if installed and declared_version != installed:
        fail(
            f"VERSION in config ({declared_version!r}) does not match the "
            f"installed MONSDA ({installed!r}). MONSDA refuses to run on a "
            "mismatch.",
            errors,
        )


def resolve(base, path):
    if not path:
        return None
    expanded = os.path.expandvars(os.path.expanduser(path))
    return expanded if os.path.isabs(expanded) else os.path.join(base, expanded)


def check_leaf_paths(leaf, base, label, errors, warnings):
    for key in ("REFERENCE", "INDEX"):
        val = leaf.get(key)
        if isinstance(val, str) and val:
            full = resolve(base, val)
            if not os.path.exists(full):
                warnings.append(f"{label}: {key} path not found: {val}")
    ann = leaf.get("ANNOTATION")
    if isinstance(ann, dict):
        for key, val in ann.items():
            if isinstance(val, str) and val:
                full = resolve(base, val)
                if not os.path.exists(full):
                    warnings.append(f"{label}: ANNOTATION.{key} path not found: {val}")
    samples = leaf.get("SAMPLES")
    if not isinstance(samples, list) or not samples:
        fail(f"{label}: SAMPLES is missing or empty", errors)
    for extra_key in ("GROUPS", "BATCHES", "TYPES"):
        extra = leaf.get(extra_key)
        if extra is not None and isinstance(samples, list) and len(extra) != len(samples):
            fail(
                f"{label}: {extra_key} has {len(extra)} entries but SAMPLES "
                f"has {len(samples)}",
                errors,
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", help="Path to the MONSDA config JSON")
    parser.add_argument(
        "--base", default=None,
        help="Directory to resolve relative paths against (default: config's own directory)",
    )
    args = parser.parse_args()

    errors, warnings = [], []

    try:
        with open(args.config) as fh:
            config = json.load(fh)
    except json.JSONDecodeError as exc:
        print(f"ERROR: {args.config} is not valid JSON: {exc}", file=sys.stderr)
        print(
            "Hint: configs/template.json has trailing '#comments' for "
            "documentation and is NOT valid JSON on its own; start from "
            "configs/template_clean.json instead.",
            file=sys.stderr,
        )
        sys.exit(2)
    except OSError as exc:
        print(f"ERROR: could not read {args.config}: {exc}", file=sys.stderr)
        sys.exit(2)

    base = args.base or os.path.dirname(os.path.abspath(args.config)) or "."

    for required in ("WORKFLOWS", "VERSION", "SETTINGS"):
        if required not in config:
            fail(f"Missing top-level key: {required}", errors)

    workflows_raw = config.get("WORKFLOWS", "")
    steps = [s.strip() for s in str(workflows_raw).split(",") if s.strip()]
    if not steps:
        fail("WORKFLOWS is empty: no workflow step selected", errors)
    for step in steps:
        if step not in KNOWN_WORKFLOW_STEPS:
            warnings.append(
                f"WORKFLOWS lists {step!r}, which is not in the known step "
                f"list ({sorted(KNOWN_WORKFLOW_STEPS)}); check reference/workflows.md"
            )
        elif step not in config:
            fail(
                f"WORKFLOWS includes {step!r} but there is no top-level "
                f"{step!r} block with its TOOLS/OPTIONS",
                errors,
            )

    version = config.get("VERSION")
    if version is not None:
        check_installed_version(version, errors, warnings)

    settings = config.get("SETTINGS")
    if isinstance(settings, dict):
        leaves = list(walk_settings_leaves(settings))
        if not leaves:
            fail("SETTINGS has no condition-tree leaf with a SAMPLES list", errors)
        for path, leaf in leaves:
            label = "SETTINGS." + ".".join(path) if path else "SETTINGS"
            check_leaf_paths(leaf, base, label, errors, warnings)

    print(f"Checked {args.config} (base dir: {base})")
    if warnings:
        print(f"\n{len(warnings)} warning(s):")
        for w in warnings:
            print(f"  - {w}")
    if errors:
        print(f"\n{len(errors)} error(s):")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    print("\nNo structural errors found. This does not replace a real "
          "`monsda --save` + dry-run (see reference/cli.md).")


if __name__ == "__main__":
    main()
