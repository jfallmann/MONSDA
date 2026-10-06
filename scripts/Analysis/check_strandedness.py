#!/usr/bin/env python3
"""RustQC strandedness runtime guard.

Compares the strand assignment inferred from a RustQC ``infer_experiment``
report against the strand setting configured in MONSDA. Runs as a gate after
the RustQC command inside the QC workflow: a confident mismatch aborts the
job (exit 2) instead of silently producing mis-stranded QC output.

Exit codes: 0 on match or when the data does not allow a confident verdict
(inconclusive/ambiguous, reported as NOT validated), 2 on a clear mismatch or
on malformed/missing input. A JSON report is written for every successfully
parsed input when ``--report`` is given.
"""

import argparse
import json
import math
import re
import sys
from pathlib import Path

PAIRED_FORWARD = "1++,1--,2+-,2-+"
PAIRED_REVERSE = "1+-,1-+,2++,2--"
SINGLE_FORWARD = "++,--"
SINGLE_REVERSE = "+-,-+"

EXPECTED_ALIASES = {
    "forward": "forward",
    "reverse": "reverse",
    "unstranded": "unstranded",
    "fr": "forward",
    "rf": "reverse",
}

THRESHOLDS = {
    "assigned_min": 0.5,
    "forward_ratio": 0.8,
    "unstranded_low": 0.4,
    "unstranded_high": 0.6,
}

SUM_TOL = 0.001


class StrandednessError(Exception):
    pass


def norm_pattern(pattern):
    return "".join(pattern.split())


def find_infer_file(qc_dir):
    matches = sorted(qc_dir.rglob("*.infer_experiment.txt"))
    if not matches:
        raise StrandednessError(
            "no *.infer_experiment.txt found under %s" % qc_dir
        )
    if len(matches) > 1:
        raise StrandednessError(
            "multiple *.infer_experiment.txt found under %s: %s"
            % (qc_dir, ", ".join(str(m) for m in matches))
        )
    return matches[0]


def parse_infer_experiment(path):
    mode = None
    fields = {}
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.strip()
            if not line:
                continue
            head = re.match(r"This is (Pair|Single)End Data$", line)
            if head:
                if mode is not None:
                    raise StrandednessError(
                        "%s:%d: duplicate paired/single header" % (path, lineno)
                    )
                mode = "paired" if head.group(1) == "Pair" else "single"
                continue
            frac = re.match(
                r'Fraction of reads (?:failed to determine|explained by "([^"]*)"):\s*(\S+)$',
                line,
            )
            if not frac:
                raise StrandednessError(
                    "%s:%d: unrecognized line %r" % (path, lineno, line)
                )
            pattern, value = frac.group(1), frac.group(2)
            key = "failed" if pattern is None else norm_pattern(pattern)
            if key not in ("failed", PAIRED_FORWARD, PAIRED_REVERSE,
                           SINGLE_FORWARD, SINGLE_REVERSE):
                raise StrandednessError(
                    "%s:%d: unknown strand pattern %r" % (path, lineno, key)
                )
            if key in fields:
                raise StrandednessError(
                    "%s:%d: duplicate field %r" % (path, lineno, key)
                )
            try:
                number = float(value)
            except ValueError:
                raise StrandednessError(
                    "%s:%d: non-numeric fraction %r" % (path, lineno, value)
                )
            if not math.isfinite(number) or not 0.0 <= number <= 1.0:
                raise StrandednessError(
                    "%s:%d: fraction %r out of [0,1]" % (path, lineno, value)
                )
            fields[key] = number
    if mode is None:
        raise StrandednessError(
            "%s: no 'This is PairEnd/SingleEnd Data' line found" % path
        )
    if "failed" not in fields:
        raise StrandednessError(
            "%s: missing 'Fraction of reads failed to determine' line" % path
        )
    if mode == "paired":
        if PAIRED_FORWARD not in fields or PAIRED_REVERSE not in fields:
            raise StrandednessError(
                "%s: paired data missing 1++/1-- pattern fractions" % path
            )
        forward = fields[PAIRED_FORWARD]
        reverse = fields[PAIRED_REVERSE]
        for key in (SINGLE_FORWARD, SINGLE_REVERSE):
            if fields.get(key, 0.0) > 0.0:
                raise StrandednessError(
                    "%s: paired data contains nonzero single-end pattern %r"
                    % (path, key)
                )
    else:
        if SINGLE_FORWARD not in fields or SINGLE_REVERSE not in fields:
            raise StrandednessError(
                "%s: single-end data missing ++/-- pattern fractions" % path
            )
        forward = fields[SINGLE_FORWARD]
        reverse = fields[SINGLE_REVERSE]
        for key in (PAIRED_FORWARD, PAIRED_REVERSE):
            if key in fields:
                raise StrandednessError(
                    "%s: single-end data contains paired pattern %r" % (path, key)
                )
    failed = fields["failed"]
    total = failed + forward + reverse
    if not (abs(total - 1.0) <= SUM_TOL or (failed == 0.0 and forward == 0.0 and reverse == 0.0)):
        raise StrandednessError(
            "%s: fractions do not sum to 1 (got %.4f)" % (path, total)
        )
    return mode, failed, forward, reverse


def infer_strandedness(failed, forward, reverse):
    assigned = forward + reverse
    if assigned < THRESHOLDS["assigned_min"]:
        return "inconclusive"
    ratio = forward / assigned
    if ratio >= THRESHOLDS["forward_ratio"]:
        return "forward"
    if reverse / assigned >= THRESHOLDS["forward_ratio"]:
        return "reverse"
    if THRESHOLDS["unstranded_low"] <= ratio <= THRESHOLDS["unstranded_high"]:
        return "unstranded"
    return "ambiguous"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Check RustQC inferred strandedness against the configured setting"
    )
    parser.add_argument("--qc-dir", required=True, help="per-sample RustQC output directory")
    parser.add_argument(
        "--expected",
        required=True,
        choices=sorted(EXPECTED_ALIASES),
        help="configured strandedness (forward/reverse/unstranded or fr/rf)",
    )
    parser.add_argument("--sample", required=True, help="sample name")
    parser.add_argument("--report", help="optional JSON report path")
    args = parser.parse_args(argv)

    expected = EXPECTED_ALIASES[args.expected]
    try:
        infer_file = find_infer_file(Path(args.qc_dir))
        mode, failed, forward, reverse = parse_infer_experiment(infer_file)
    except StrandednessError as err:
        print("ERROR: %s" % err, file=sys.stderr)
        return 2

    inferred = infer_strandedness(failed, forward, reverse)
    if inferred == expected:
        status = "match"
    elif inferred in ("forward", "reverse", "unstranded"):
        status = "mismatch"
    else:
        status = inferred

    report = {
        "sample": args.sample,
        "expected": expected,
        "inferred": inferred,
        "status": status,
        "paired": mode == "paired",
        "fractions": {
            "failed": failed,
            "forward": forward,
            "reverse": reverse,
        },
        "thresholds": dict(THRESHOLDS),
    }
    if args.report:
        Path(args.report).write_text(json.dumps(report, indent=2) + "\n")

    if status == "match":
        return 0
    if status == "mismatch":
        print(
            "ERROR: strandedness mismatch for sample %s: configured %s, inferred %s "
            "(forward=%.4f, reverse=%.4f, failed=%.4f). Correct "
            "SETTINGS...SEQUENCING (paired,fr/rf or single,fr/rf) and check "
            "protocol/reference. Not auto-changing."
            % (args.sample, expected, inferred, forward, reverse, failed),
            file=sys.stderr,
        )
        return 2
    print(
        "WARNING: strandedness NOT validated for sample %s: expected %s, inferred %s "
        "(forward=%.4f, reverse=%.4f, failed=%.4f)"
        % (args.sample, expected, inferred, forward, reverse, failed),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
