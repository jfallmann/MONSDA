#! /usr/bin/env python3

"""
Normalize an oarfish `.quant` file (tab-separated columns:
tname, len, num_reads - no header fields for TPM/EffectiveLength) into a
salmon-compatible `quant.sf` so downstream DTU tooling (tximport type="salmon",
build_DTU_table.py, SPIT.py) can consume long-read quantification output
without any further changes.

oarfish does not report an effective length (no fragment-length bias model
for long reads), so transcript length is used as the effective length, as
recommended by the oarfish authors for downstream tximport usage.
"""

import argparse
import csv
import gzip
import inspect
import logging
import os
import sys
import traceback as tb

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

scriptname = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")
log = logging.getLogger(scriptname)


def parseargs():
    parser = argparse.ArgumentParser(
        description="Convert an oarfish .quant file to a salmon-compatible quant.sf"
    )
    parser.add_argument(
        "-i", "--input", required=True, type=str, help="oarfish .quant input file"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=str,
        help="Output path for the salmon-style quant.sf(.gz)",
    )
    parser.add_argument("--loglevel", default="INFO", help="Log verbosity")
    return parser.parse_args()


def convert(inpath, outpath):
    logid = scriptname + ".convert: "

    rows = []
    with open(inpath, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = {name.lower(): name for name in reader.fieldnames or []}
        for required in ("tname", "len", "num_reads"):
            if required not in fieldnames:
                raise ValueError(
                    f"{inpath} is missing required oarfish column '{required}', "
                    f"found columns: {reader.fieldnames}"
                )
        for row in reader:
            name = row[fieldnames["tname"]]
            length = float(row[fieldnames["len"]])
            num_reads = float(row[fieldnames["num_reads"]])
            rows.append((name, length, num_reads))

    # TPM: standard two-pass rate normalization, using transcript length as
    # the effective length since oarfish reports none.
    rates = [
        (num_reads / length if length > 0 else 0.0) for _, length, num_reads in rows
    ]
    rate_sum = sum(rates)
    scale = 1e6 / rate_sum if rate_sum > 0 else 0.0

    opn = gzip.open if str(outpath).endswith(".gz") else open
    with opn(outpath, "wt") as out:
        out.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        for (name, length, num_reads), rate in zip(rows, rates):
            tpm = rate * scale
            out.write(
                f"{name}\t{int(length)}\t{length:.3f}\t{tpm:.6f}\t{num_reads:.3f}\n"
            )

    log.info(logid + f"Wrote {len(rows)} transcripts from {inpath} to {outpath}")


if __name__ == "__main__":
    logid = scriptname + ".main: "
    try:
        args = parseargs()
        try:
            log = setup_logger(
                name=scriptname,
                log_file="stderr",
                logformat="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
                datefmt="%m-%d %H:%M",
                level=args.loglevel,
            )
        except Exception:
            log = logging.getLogger(scriptname)
        convert(args.input, args.output)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(exc_type, exc_value, exc_tb)
        log.error(logid + "".join(tbe.format()))
        sys.exit(1)

#
# oarfish_to_salmon.py ends here
