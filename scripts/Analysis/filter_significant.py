#!/usr/bin/env python3
"""Header-aware significant-hit filter for DE/DEU/DAS/DTU result tables.

Replaces the positional Perl one-liners in the workflows. Columns are located
by header name, so adding or reordering annotation columns (Coordinates, Gene,
stat, ...) cannot silently shift the tested fields. Different tools name their
columns differently (log2FoldChange/logFC/lfc, padj/FDR/adj_pvalue), so the
caller passes the names of the table it filters.

Numeric mode (--effect-column/--adjusted-p-column): a row is kept when the
adjusted p-value is strictly below --p-cutoff and the absolute effect size is
at least --lfc-cutoff. Flag mode (--flag-column/--flag-value): a row is kept
when the flag column equals the given value, used for tools such as DIEGO that
report significance categorically instead of as p-value plus effect size. Both
modes can be combined, a row then has to satisfy all criteria.

Only finite values in the tested numeric columns are evaluated; missing
optional annotation columns or an NA in an untested column do not discard a
row. Non-numeric values in a tested column are treated as malformed input and
abort with a clear error. Header-only inputs produce outputs containing only
the header. Outputs are gzip compressed when their name ends in .gz, plain text
otherwise.
"""

import argparse
import contextlib
import gzip
import math
import sys

NON_FINITE = {"", "NA", "N/A", "NaN", "nan", "Inf", "-Inf", "inf", "-inf"}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def open_out(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", encoding="utf-8", newline="")
    return open(path, "w", encoding="utf-8", newline="")


def column_index(header, colname, what):
    if colname not in header:
        raise ValueError(
            "%s column %r not found in header: %s" % (what, colname, header)
        )
    return header.index(colname)


def read_header(handle):
    line = handle.readline()
    if not line:
        raise ValueError("input file is empty, no header line found")
    header = line.rstrip("\n").split("\t")
    if not header or any(not col for col in header):
        raise ValueError("malformed header line")
    return header


def parse_required(value, colname, lineno):
    stripped = value.strip()
    if stripped in NON_FINITE:
        return None
    try:
        return float(stripped)
    except ValueError:
        raise ValueError(
            "line %d: non-numeric value %r in required column %r"
            % (lineno, value, colname)
        )


def filter_table(inpath, effect_col, adjp_col, p_cut, lfc_cut, flag_col, flag_value,
                 out_sig, out_up, out_down):
    with open_maybe_gzip(inpath) as fh, contextlib.ExitStack() as stack:
        header = read_header(fh)
        eidx = aidx = fidx = None
        if effect_col is not None:
            eidx = column_index(header, effect_col, "effect")
            aidx = column_index(header, adjp_col, "adjusted p")
        if flag_col is not None:
            fidx = column_index(header, flag_col, "flag")
        header_line = "\t".join(header) + "\n"
        sig = stack.enter_context(open_out(out_sig))
        up = stack.enter_context(open_out(out_up)) if out_up else None
        down = stack.enter_context(open_out(out_down)) if out_down else None
        for handle in (sig, up, down):
            if handle is not None:
                handle.write(header_line)
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(
                    "line %d: expected %d fields, got %d"
                    % (lineno, len(header), len(fields))
                )
            if fidx is not None and fields[fidx].strip() != flag_value:
                continue
            lfc = None
            if eidx is not None:
                pval = parse_required(fields[aidx], adjp_col, lineno)
                lfc = parse_required(fields[eidx], effect_col, lineno)
                if pval is None or lfc is None:
                    continue
                if not (math.isfinite(pval) and math.isfinite(lfc)):
                    continue
                if not (pval < p_cut and abs(lfc) >= lfc_cut):
                    continue
            sig.write(line)
            if lfc is not None:
                if up is not None and lfc >= lfc_cut:
                    up.write(line)
                if down is not None and lfc <= -lfc_cut:
                    down.write(line)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter a DE result table for significant hits by column name."
    )
    parser.add_argument("--effect-column", help="effect size column (log2FoldChange/logFC/lfc)")
    parser.add_argument("--adjusted-p-column", help="adjusted p-value column (padj/FDR/adj_pvalue)")
    parser.add_argument("--p-cutoff", type=float, help="adjusted p-value cutoff (strict <)")
    parser.add_argument("--lfc-cutoff", type=float, help="absolute log fold change cutoff (>=)")
    parser.add_argument("--flag-column", help="categorical significance column (DIEGO: significant)")
    parser.add_argument("--flag-value", default="yes", help="value the flag column has to equal")
    parser.add_argument("--input", required=True, help="input result table (plain or gzip)")
    parser.add_argument("--output-sig", required=True, help="output with all significant rows")
    parser.add_argument("--output-up", help="output with up-regulated rows, numeric mode only")
    parser.add_argument("--output-down", help="output with down-regulated rows, numeric mode only")
    args = parser.parse_args(argv)
    numeric = [args.effect_column, args.adjusted_p_column, args.p_cutoff, args.lfc_cutoff]
    if any(v is not None for v in numeric) and any(v is None for v in numeric):
        parser.error(
            "--effect-column, --adjusted-p-column, --p-cutoff and --lfc-cutoff "
            "have to be given together"
        )
    if args.effect_column is None and args.flag_column is None:
        parser.error("either the numeric columns or --flag-column has to be given")
    if args.effect_column is None and (args.output_up or args.output_down):
        parser.error("--output-up/--output-down need --effect-column")
    try:
        filter_table(
            args.input,
            args.effect_column,
            args.adjusted_p_column,
            args.p_cutoff,
            args.lfc_cutoff,
            args.flag_column,
            args.flag_value,
            args.output_sig,
            args.output_up,
            args.output_down,
        )
    except (OSError, ValueError) as exc:
        print("filter_significant: %s" % exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
