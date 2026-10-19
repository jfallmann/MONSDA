#! /usr/bin/env python3
"""Aggregate oarfish transcript-level long-read quantification to gene level
and emit a featureCounts-compatible count file, so the existing
DE build_count_table.py/DESEQ2.R/EDGER.R pipeline can consume it unchanged.

oarfish's ``.quant`` output has no TPM/EffectiveLength, just:
    tname\tlen\tnum_reads
Gene coordinates and the transcript->gene mapping are taken from the GTF/GFF
annotation used for mapping/quantification.
"""

import sys
import os
import re
import gzip
import argparse
import inspect

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *
import logging

log = logging.getLogger(os.path.basename(__file__).replace(".py", ""))


def parseargs():
    parser = argparse.ArgumentParser(
        description="Aggregate oarfish transcript counts to gene level in featureCounts format"
    )
    parser.add_argument(
        "-i", "--quant", required=True, help="oarfish .quant output file"
    )
    parser.add_argument(
        "-g", "--gtf", required=True, help="Annotation GTF/GFF used for quantification"
    )
    parser.add_argument(
        "-s",
        "--sample",
        required=True,
        help="Sample column name to use in the output header (e.g. the bam filename)",
    )
    parser.add_argument("-o", "--output", required=True, help="Output count file (.gz)")
    parser.add_argument("--loglevel", default="INFO", help="Log verbosity")
    return parser.parse_args()


def openfile(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_gtf_attrs(field):
    attrs = {}
    if "=" in field and '"' not in field:
        for kv in field.strip().split(";"):
            kv = kv.strip()
            if not kv or "=" not in kv:
                continue
            k, v = kv.split("=", 1)
            attrs[k.strip()] = v.strip()
    else:
        for m in _ATTR_RE.finditer(field):
            attrs[m.group(1)] = m.group(2)
    return attrs


def build_tx2gene_and_coords(gtf):
    """Parse transcript->gene mapping and gene coordinates (min start, max
    end across all exon/transcript lines seen per gene) from a GTF/GFF."""
    tx2gene = {}
    gene_coords = {}
    with openfile(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrfield = cols[:9]
            attrs = parse_gtf_attrs(attrfield)
            gene_id = attrs.get("gene_id") or attrs.get("gene")
            tx_id = attrs.get("transcript_id")
            if tx_id and gene_id and tx_id not in tx2gene:
                tx2gene[tx_id] = gene_id
            if gene_id:
                s, e = int(start), int(end)
                if gene_id in gene_coords:
                    oc, os_, oe, ost = gene_coords[gene_id]
                    gene_coords[gene_id] = (oc, min(os_, s), max(oe, e), ost)
                else:
                    gene_coords[gene_id] = (chrom, s, e, strand)
    return tx2gene, gene_coords


def aggregate(quant, tx2gene):
    gene_counts = {}
    unmapped = 0
    with openfile(quant) as fh:
        header = fh.readline()
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            tname, _tlen, num_reads = fields[0], fields[1], fields[2]
            gene_id = tx2gene.get(tname)
            if not gene_id:
                unmapped += 1
                continue
            gene_counts[gene_id] = gene_counts.get(gene_id, 0.0) + float(num_reads)
    if unmapped:
        log.warning(
            f"oarfish_to_featurecounts: {unmapped} transcript(s) from {quant} had no gene_id in the annotation and were skipped"
        )
    return gene_counts


def write_featurecounts(output, sample, gene_counts, gene_coords):
    with gzip.open(output, "wt") as out:
        out.write(f"# Program:oarfish; Command:aggregated to gene level by oarfish_to_featurecounts.py\n")
        out.write("Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + sample + "\n")
        for gene_id in sorted(gene_counts):
            chrom, start, end, strand = gene_coords.get(gene_id, ("NA", 0, 0, "."))
            length = max(end - start + 1, 1)
            out.write(
                f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\t{length}\t{round(gene_counts[gene_id])}\n"
            )


def main():
    args = parseargs()
    log.setLevel(args.loglevel)
    tx2gene, gene_coords = build_tx2gene_and_coords(args.gtf)
    gene_counts = aggregate(args.quant, tx2gene)
    write_featurecounts(args.output, args.sample, gene_counts, gene_coords)
    log.info(f"oarfish_to_featurecounts: wrote {len(gene_counts)} genes to {args.output}")


if __name__ == "__main__":
    main()
