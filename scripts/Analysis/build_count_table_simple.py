#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import sys
import traceback as tb


def parseargs():
    parser = argparse.ArgumentParser(
        description="Simple script to build a count table from replicate files"
    )
    parser.add_argument(
        "-r",
        "--replicates",
        required=True,
        type=str,
        help="Comma-separated list of replicate count files",
    )
    parser.add_argument(
        "--table",
        dest="table",
        required=True,
        type=str,
        default="counts.table",
        help="Name of table to write to",
    )
    parser.add_argument(
        "--anno",
        dest="anno",
        required=True,
        type=str,
        default="counts.anno",
        help="Name of anno to write to",
    )
    parser.add_argument("--loglevel", default="INFO", help="Log verbosity")
    return parser.parse_args()


def setup_logger(loglevel):
    log = logging.getLogger("build_count_table_simple")
    lvl = getattr(logging, loglevel.upper(), logging.INFO)
    log.setLevel(lvl)
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(message)s", datefmt="%m-%d %H:%M"
        )
    )
    log.addHandler(handler)
    return log


def prepare_table(replicates, table, anno, log):
    try:
        replist = [r for r in replicates.strip().split(",") if r]
        log.info(f"Replicates: {replist}")

        myMatrix = []
        myMatrix.append(["ID"])  # header row
        myMatrix[0].append("Length")
        for idx, rep in enumerate(replist):
            rep_name = os.path.basename(rep).split("_mapped")[0]
            myMatrix[0].append(rep_name)
            if ".gz" in rep:
                myInput = gzip.open(rep, "rt")
            else:
                myInput = open(rep, "r")

            lineNumber = 0
            for line in myInput:
                if "#" in line[0:5] or ".bam" in line[-10:]:
                    continue
                columns = line.strip().split("\t")
                if (
                    columns[0] != "name"
                    and columns[0] != "Geneid"
                    and columns[1] != "count"
                ):
                    lineNumber += 1
                    if idx == 0:
                        myMatrix.append([columns[0], columns[5]])
                    if len(columns) > 1 and columns[-1] != columns[0]:
                        myMatrix[lineNumber].append(round(float(columns[-1])))
                    else:
                        myMatrix[lineNumber].append("0")
            myInput.close()

        # Write annotation file (just sample names)
        annos = [name for name in myMatrix[0][1:]]
        with gzip.open(anno, "wt") as a:
            for c in annos:
                a.write(f"{c}\n") if c != "Length" else None

        # Write count table
        with gzip.open(table, "wt") as t:
            for row in myMatrix:
                t.write("\t".join(map(str, row)) + "\n")

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(exc_type, exc_value, exc_tb)
        log.error("".join(tbe.format()))


if __name__ == "__main__":
    args = parseargs()
    log = setup_logger(args.loglevel)
    prepare_table(args.replicates, args.table, args.anno, log)
