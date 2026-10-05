#!/usr/bin/env python3
"""Header-aware significant-hit filter for DE result tables.

Replaces the positional Perl one-liners in the deseq2_DE/edger_DE workflows.
Columns are located by header name, so adding or reordering annotation columns
(Coordinates, stat, ...) cannot silently shift the tested fields.

A row is kept when the adjusted p-value is strictly below --p-cutoff and the
absolute effect size is at least --lfc-cutoff. Only finite values in the two
required columns are evaluated; missing optional annotation columns (e.g.
Coordinates, Gene) or an NA in an optional column such as edgeR's F do not
discard a row. Non-numeric values in a required column are treated as malformed
input and abort with a clear error. Header-only inputs produce gzip outputs
containing only the header.
"""

import argparse
import gzip
import math
import sys

NON_FINITE = {"", "NA", "N/A", "NaN", "nan", "Inf", "-Inf", "inf", "-inf"}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


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


def filter_table(inpath, effect_col, adjp_col, p_cut, lfc_cut, out_sig, out_up, out_down):
    with open_maybe_gzip(inpath) as fh:
        header = read_header(fh)
        if effect_col not in header:
            raise ValueError(
                "effect column %r not found in header: %s" % (effect_col, header)
            )
        if adjp_col not in header:
            raise ValueError(
                "adjusted p column %r not found in header: %s" % (adjp_col, header)
            )
        eidx = header.index(effect_col)
        aidx = header.index(adjp_col)
        header_line = "\t".join(header) + "\n"
        with gzip.open(out_sig, "wt", encoding="utf-8", newline="") as sig, \
             gzip.open(out_up, "wt", encoding="utf-8", newline="") as up, \
             gzip.open(out_down, "wt", encoding="utf-8", newline="") as down:
            sig.write(header_line)
            up.write(header_line)
            down.write(header_line)
            for lineno, line in enumerate(fh, start=2):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != len(header):
                    raise ValueError(
                        "line %d: expected %d fields, got %d"
                        % (lineno, len(header), len(fields))
                    )
                pval = parse_required(fields[aidx], adjp_col, lineno)
                lfc = parse_required(fields[eidx], effect_col, lineno)
                if pval is None or lfc is None:
                    continue
                if not (math.isfinite(pval) and math.isfinite(lfc)):
                    continue
                if pval < p_cut and abs(lfc) >= lfc_cut:
                    sig.write(line)
                    if lfc >= lfc_cut:
                        up.write(line)
                    if lfc <= -lfc_cut:
                        down.write(line)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter a DE result table for significant hits by column name."
    )
    parser.add_argument("--effect-column", required=True, help="effect size column (log2FoldChange/logFC)")
    parser.add_argument("--adjusted-p-column", required=True, help="adjusted p-value column (padj/FDR)")
    parser.add_argument("--p-cutoff", required=True, type=float, help="adjusted p-value cutoff (strict <)")
    parser.add_argument("--lfc-cutoff", required=True, type=float, help="absolute log fold change cutoff (>=)")
    parser.add_argument("--input", required=True, help="input result table (plain or gzip)")
    parser.add_argument("--output-sig", required=True, help="gzip output with all significant rows")
    parser.add_argument("--output-up", required=True, help="gzip output with up-regulated rows")
    parser.add_argument("--output-down", required=True, help="gzip output with down-regulated rows")
    args = parser.parse_args(argv)
    try:
        filter_table(
            args.input,
            args.effect_column,
            args.adjusted_p_column,
            args.p_cutoff,
            args.lfc_cutoff,
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
