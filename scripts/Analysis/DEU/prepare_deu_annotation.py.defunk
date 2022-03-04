#!/usr/bin/env python3

import glob
import os
import sys
import inspect
import json
import shutil
from collections import defaultdict
import traceback as tb
import argparse
import re
import gzip
import logging

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")
    log = logging.getLogger(scriptname)

    lvl = log.level if log.level else "INFO"
    for handler in log.handlers[:]:
        handler.close()
        log.removeHandler(handler)

    handler = logging.FileHandler("LOGS/MONSDA.log", mode="a")
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)-8s %(name)-12s %(message)s",
            datefmt="%m-%d %H:%M",
        )
    )
    log.addHandler(handler)
    log.setLevel(lvl)

except Exception:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type,
        exc_value,
        exc_tb,
    )
    print("".join(tbe.format()), file=sys.stderr)


def parseargs():
    parser = argparse.ArgumentParser(
        description="Converting ENSEMBL gff3 to DEXSeq ready gff3"
    )
    parser.add_argument(
        "-g", "--gff", required=True, help="gff3 input file in ENSEMBL format"
    )
    parser.add_argument("-o", "--outgff", required=True, help="Output gff3 to write to")
    parser.add_argument(
        "-v",
        "--loglevel",
        type=str,
        default="INFO",
        choices=["WARNING", "ERROR", "INFO", "DEBUG"],
        help="Set log level",
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()


def reformat(gff=None, outgff=None):
    try:
        print(gff, outgff)
        if os.path.isfile(os.path.abspath(gff)):
            if ".gz" in gff:
                f = gzip.open(gff, "rt")
            else:
                f = open(gff, "r")
        else:
            f = gff

        mapping = defaultdict()
        out = list()

        for line in f:
            line = line.rstrip()
            if "#" in line[:3]:
                out.append(line)
            else:
                entries = line.rstrip().split("\t")
                anno = entries[-1]
                ids = anno.split(";")
                if entries[2] != "exon" and "Parent" in anno:
                    trans, parent = None, None
                    for id in ids:
                        if "transcript:" in id:
                            trans = id.split(":")[-1]
                        elif "gene:" in id:
                            parent = id.split(":")[-1]

                    if trans and parent:
                        mapping[trans] = parent
                    out.append(line)

                elif entries[2] == "exon":
                    parenttrans = ids[0].split(":")[-1]
                    line += ";ParentGene=" + mapping[parenttrans]
                    out.append(line)

                else:
                    out.append(line)

        if ".gz" in outgff:
            with gzip.open(outgff, "wb") as output:
                output.write(bytes(str.join("\n", out), encoding="UTF-8"))
        else:
            with open(outgff, "w") as output:
                output.write(str.join("\n", out))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )


####################
####    MAIN    ####
####################

if __name__ == "__main__":

    try:
        args = parseargs()
        knownargs = args[0]
        reformat(knownargs.gff, knownargs.outgff)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
