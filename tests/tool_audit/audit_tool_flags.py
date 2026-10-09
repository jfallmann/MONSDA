#!/usr/bin/env python3
"""Audit MONSDA's hardcoded external-tool CLI flags against reality.

Runs each tool named in tool_manifest.yaml inside its apptainer container
and checks that the flags MONSDA hardcodes in workflows/*.smk and
workflows/*.nf are still accepted by the tool's own --help (or, for tools
whose --help exits before printing usage, a minimal real run).

Usage:
    conda activate apptainer
    python tests/tool_audit/audit_tool_flags.py [--strict] [--tool NAME ...]

Requires:
    - apptainer on PATH (e.g. `conda activate apptainer`)
    - container .sif files under --container-dir (default matches the
      user's local ~/Work2/Container/MONSDA)

See README.md in this directory for the full manual findings this manifest
was derived from.
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.stderr.write(
        "error: PyYAML is required (pip install pyyaml / conda install pyyaml)\n"
    )
    sys.exit(2)

HERE = Path(__file__).resolve().parent
DEFAULT_MANIFEST = HERE / "tool_manifest.yaml"
DEFAULT_CONTAINER_DIR = Path.home() / "Work2" / "Container" / "MONSDA"
DEFAULT_CONTAINER_SUFFIX = "-1.5.0.sif"
DEFAULT_TIMEOUT = 60

STATUS_PASS = "PASS"
STATUS_DRIFT = "DRIFT"
STATUS_GAP = "GAP"
STATUS_SKIP = "SKIP"


class Result:
    def __init__(self, tool, invocation, status, detail):
        self.tool = tool
        self.invocation = invocation
        self.status = status
        self.detail = detail

    def as_dict(self):
        return {
            "tool": self.tool,
            "invocation": self.invocation,
            "status": self.status,
            "detail": self.detail,
        }


def load_manifest(path):
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return data.get("tools", {})


def apptainer_available():
    return shutil.which("apptainer") is not None


def container_path(container_dir, container_name, suffix):
    return Path(container_dir) / f"{container_name}{suffix}"


def run_apptainer(sif, binary, args, timeout, extra_binds=None):
    cmd = ["apptainer", "exec"]
    for b in extra_binds or []:
        cmd += ["--bind", b]
    cmd += [str(sif), binary, *args]
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout
        )
        return proc.returncode, proc.stdout, proc.stderr, None
    except subprocess.TimeoutExpired:
        return None, "", "", f"timed out after {timeout}s"
    except FileNotFoundError as e:
        return None, "", "", str(e)


def binary_missing(returncode, stderr):
    if returncode == 127:
        return True
    low = stderr.lower()
    return "executable file not found" in low or "command not found" in low


def check_help_invocation(tool, container_name, binary, sif, inv, timeout):
    name = inv.get("name", "help")
    args = inv.get("args", ["--help"])
    expect_exit = inv.get("expect_exit", [0])
    expected_flags = inv.get("expected_flags", [])
    expected_text = inv.get("expected_text", [])
    expect_value_required = inv.get("expect_value_required", False)

    rc, out, err, launch_err = run_apptainer(sif, binary, args, timeout)
    if launch_err is not None:
        return Result(tool, name, STATUS_GAP, f"could not launch: {launch_err}")
    if binary_missing(rc, err):
        return Result(
            tool, name, STATUS_GAP, f"binary '{binary}' not found in container"
        )

    combined = out + "\n" + err
    if rc not in expect_exit:
        return Result(
            tool,
            name,
            STATUS_DRIFT,
            f"exit code {rc} not in expected {expect_exit}; output head: "
            f"{combined.strip()[:200]!r}",
        )

    missing = [f for f in expected_flags if f not in combined]
    missing += [t for t in expected_text if t not in combined]
    if missing:
        return Result(
            tool, name, STATUS_DRIFT, f"flags/text missing from help output: {missing}"
        )

    if expect_value_required:
        # Regression check: the flag (first entry in expected_flags) must
        # error out if invoked bare/without a value, not silently accept it.
        flag = expected_flags[0]
        rc2, out2, err2, launch_err2 = run_apptainer(
            sif, binary, args[:-1] + [flag] if args and args[-1] == "--help" else [flag],
            timeout,
        )
        # Fall back: just invoke the subcommand (without --help) plus the bare flag
        if launch_err2 is None and not binary_missing(rc2, err2):
            combined2 = (out2 + err2).lower()
            if rc2 == 0 or "value is required" not in combined2 and "requires a value" not in combined2 and "requires an argument" not in combined2:
                return Result(
                    tool,
                    name,
                    STATUS_DRIFT,
                    f"expected '{flag}' to require an explicit value again "
                    "(regression check for old bare-flag behaviour), but "
                    f"invocation did not error as expected (exit={rc2})",
                )

    return Result(tool, name, STATUS_PASS, "ok")


def run_prep_steps(prep, tmpdir, fixture_path):
    """Execute manifest-declared `prep` steps (write/gzip) before a
    minimal_run invocation, so fixtures needing more than one static file
    (e.g. the oarfish gzipped-GTF regression probe) don't need to be
    checked into fixtures/ as binary blobs."""
    for step in prep or []:
        if "write" in step:
            target = Path(step["write"].format(fixture=str(fixture_path) if fixture_path else "", tmpdir=tmpdir))
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(step.get("content", ""))
        elif "gzip" in step:
            import gzip as gzip_module

            src = Path(step["gzip"].format(fixture=str(fixture_path) if fixture_path else "", tmpdir=tmpdir))
            dest = Path(step["dest"].format(fixture=str(fixture_path) if fixture_path else "", tmpdir=tmpdir))
            with open(src, "rb") as fh_in, gzip_module.open(dest, "wb") as fh_out:
                fh_out.write(fh_in.read())
        else:
            raise ValueError(f"unknown prep step: {step}")


def check_minimal_run_invocation(tool, container_name, binary, sif, inv, timeout, fixtures_dir):
    name = inv.get("name", "minimal_run")
    fixture = inv.get("fixture")
    expect_exit = inv.get("expect_exit", [0])
    args_template = inv.get("args_template", [])
    expect_output_exists = inv.get("expect_output_exists")
    expect_stderr_contains = inv.get("expect_stderr_contains", [])
    prep = inv.get("prep")

    fixture_path = fixtures_dir / fixture if fixture else None
    if fixture_path is not None and not fixture_path.exists():
        return Result(
            tool, name, STATUS_GAP, f"fixture missing: {fixture_path}"
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            run_prep_steps(prep, tmpdir, fixture_path)
        except Exception as e:
            return Result(tool, name, STATUS_GAP, f"prep step failed: {e}")

        args = [
            a.format(fixture=str(fixture_path) if fixture_path else "", tmpdir=tmpdir)
            for a in args_template
        ]
        binds = [f"{tmpdir}:{tmpdir}"]
        if fixture_path is not None:
            binds.append(f"{fixture_path.parent}:{fixture_path.parent}")

        rc, out, err, launch_err = run_apptainer(sif, binary, args, timeout, extra_binds=binds)
        if launch_err is not None:
            return Result(tool, name, STATUS_GAP, f"could not launch: {launch_err}")
        if binary_missing(rc, err):
            return Result(
                tool, name, STATUS_GAP, f"binary '{binary}' not found in container"
            )
        if rc not in expect_exit:
            return Result(
                tool,
                name,
                STATUS_DRIFT,
                f"exit code {rc} not in expected {expect_exit}; stderr head: "
                f"{err.strip()[:200]!r}",
            )
        if expect_output_exists:
            out_path = Path(expect_output_exists.format(fixture=str(fixture_path) if fixture_path else "", tmpdir=tmpdir))
            if not out_path.exists() or out_path.stat().st_size == 0:
                return Result(
                    tool, name, STATUS_DRIFT, f"expected output not created: {out_path}"
                )

        if expect_stderr_contains:
            missing_text = [t for t in expect_stderr_contains if t not in err]
            if missing_text:
                return Result(
                    tool,
                    name,
                    STATUS_DRIFT,
                    f"expected stderr text not found (message may have changed "
                    f"or the underlying issue may be fixed upstream): {missing_text}; "
                    f"stderr head: {err.strip()[:300]!r}",
                )

    return Result(tool, name, STATUS_PASS, "ok")


def audit_tool(tool, spec, container_dir, suffix, timeout, fixtures_dir):
    if spec.get("na"):
        return [Result(tool, "n/a", STATUS_SKIP, spec.get("na_reason", "no hardcoded flags"))]

    container_name = spec.get("container", tool)
    binary = spec.get("binary", tool)
    sif = container_path(container_dir, container_name, suffix)

    known_gap = spec.get("known_gap", False)

    if not sif.exists():
        detail = f"container file missing: {sif}"
        return [Result(tool, "n/a", STATUS_GAP, detail)]

    results = []
    invocations = spec.get("invocations", [])
    if not invocations and known_gap:
        # Probe the binary directly to confirm/refresh the known gap.
        rc, out, err, launch_err = run_apptainer(sif, binary, ["--help"], timeout)
        if launch_err is not None or binary_missing(rc, err):
            results.append(
                Result(tool, "presence", STATUS_GAP, spec.get("gap_reason", "known gap"))
            )
        else:
            results.append(
                Result(
                    tool,
                    "presence",
                    STATUS_DRIFT,
                    "binary now FOUND in container but manifest still marks "
                    "known_gap: true and has no invocations to check - "
                    "update tool_manifest.yaml",
                )
            )
        return results

    for inv in invocations:
        if inv.get("skip"):
            results.append(
                Result(
                    tool,
                    inv.get("name", "skip"),
                    STATUS_SKIP,
                    inv.get("skip_reason", "marked skip in manifest"),
                )
            )
            continue
        mode = inv.get("mode", "help")
        if mode == "minimal_run":
            results.append(
                check_minimal_run_invocation(
                    tool, container_name, binary, sif, inv, timeout, fixtures_dir
                )
            )
        else:
            results.append(
                check_help_invocation(tool, container_name, binary, sif, inv, timeout)
            )

    # Reclassify GAP results for tools already known to be gap-y, unless
    # they newly started passing (handled above) or the gap moved (still a
    # GAP either way, just flagged non-fatal below in the summary/exit code
    # logic via `known_gap`).
    return results


def print_table(all_results, manifest):
    widths = (14, 22, 8)
    header = f"{'TOOL':<{widths[0]}} {'INVOCATION':<{widths[1]}} {'STATUS':<{widths[2]}} DETAIL"
    print(header)
    print("-" * len(header))
    for r in all_results:
        print(f"{r.tool:<{widths[0]}} {r.invocation:<{widths[1]}} {r.status:<{widths[2]}} {r.detail}")


def summarize(all_results, manifest, strict):
    counts = {STATUS_PASS: 0, STATUS_DRIFT: 0, STATUS_GAP: 0, STATUS_SKIP: 0}
    fatal = []
    for r in all_results:
        counts[r.status] = counts.get(r.status, 0) + 1
        if r.status == STATUS_DRIFT:
            fatal.append(r)
        elif r.status == STATUS_GAP:
            known = manifest.get(r.tool, {}).get("known_gap", False)
            if not known:
                fatal.append(r)
            elif strict:
                fatal.append(r)

    print()
    print(
        f"Summary: {counts[STATUS_PASS]} passed, {counts[STATUS_DRIFT]} drifted, "
        f"{counts[STATUS_GAP]} gaps, {counts[STATUS_SKIP]} skipped"
    )
    if fatal:
        print(f"FAILED: {len(fatal)} finding(s) require attention:")
        for r in fatal:
            print(f"  - [{r.status}] {r.tool} / {r.invocation}: {r.detail}")
    return fatal


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--container-dir", default=str(DEFAULT_CONTAINER_DIR))
    ap.add_argument("--container-suffix", default=DEFAULT_CONTAINER_SUFFIX)
    ap.add_argument("--tool", action="append", dest="tools", default=None,
                     help="limit to specific tool name(s); repeatable")
    ap.add_argument("--strict", action="store_true",
                     help="treat known gaps as failures too")
    ap.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT)
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)

    if not apptainer_available():
        sys.stderr.write(
            "error: 'apptainer' not found on PATH. Activate the apptainer "
            "conda env first, e.g.:\n"
            "  source ~/anaconda3/etc/profile.d/conda.sh && conda activate apptainer\n"
        )
        return 3

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        sys.stderr.write(f"error: manifest not found: {manifest_path}\n")
        return 3
    manifest = load_manifest(manifest_path)

    container_dir = Path(args.container_dir).expanduser()
    if not container_dir.exists():
        sys.stderr.write(
            f"error: container dir not found: {container_dir} "
            "(pass --container-dir to override)\n"
        )
        return 3

    fixtures_dir = HERE / "fixtures"

    tool_names = args.tools or sorted(manifest.keys())
    unknown = [t for t in tool_names if t not in manifest]
    if unknown:
        sys.stderr.write(f"error: unknown tool(s) in --tool: {unknown}\n")
        return 3

    all_results = []
    for tool in tool_names:
        spec = manifest[tool]
        all_results.extend(
            audit_tool(tool, spec, container_dir, args.container_suffix, args.timeout, fixtures_dir)
        )

    print_table(all_results, manifest)
    fatal = summarize(all_results, manifest, args.strict)

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump([r.as_dict() for r in all_results], fh, indent=2)

    return 1 if fatal else 0


if __name__ == "__main__":
    sys.exit(main())
