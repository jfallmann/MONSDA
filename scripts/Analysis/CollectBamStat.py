# CollectBamStat.py ---
#
# Filename: CollectBamStat.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue May 15 12:44:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Jun 20 11:44:19 2018 (+0200)
#           By: Joerg Fallmann
#     Update #: 320
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
#
#
#
#

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Code:

### IMPORTS
import argparse
import pprint  # pp = pprint.PrettyPrinter(indent=4) #pp.pprint(stuff)
from io import StringIO
import os
import gzip
import re
import sys
import pysam
from pyfaidx import Fasta

# import importlib
import multiprocessing
from multiprocessing import Manager
import traceback as tb
import inspect

# Container
import collections

# logging
import logging

####load own modules
cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

### MAIN

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
        description="Read BAM file read by read and collect mapping statistics"
    )
    parser.add_argument(
        "-f", "--fasta", type=str, help="Reference genome FASTA, needs to be indexed"
    )
    parser.add_argument(
        "-b",
        "--bams",
        type=str,
        help="Mapped reads BAM(s) comma separated, need to besorted and indexed",
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=".", help="Output directory"
    )
    parser.add_argument(
        "-z",
        "--procs",
        type=int,
        default=1,
        help="Number of parallel processed to run this job with",
    )
    parser.add_argument(
        "-v",
        "--verbosity",
        type=int,
        default=0,
        choices=[0, 1],
        help="increase output verbosity",
    )

    return parser.parse_args()


def collectstats(fasta, bams, outdir, procs, verbosity):

    if outdir:
        print("Checking or creating outdir " + str(outdir))
        if not os.path.isabs(outdir):
            outdir = os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        outdir = os.path.abspath(os.getcwd())

    ###Multithreading
    num_processes = procs or 1
    pool = multiprocessing.Pool(processes=num_processes)
    # 	result_list = []

    for bam in bams.split(","):
        fn = os.path.basename(bam)
        out = os.path.join(outdir, fn + ".stats.gz")

        if os.path.isfile(out):
            print("File " + out + " exists, will be deleted")
            os.remove(out)
        ####Check if index files available
        # 		try:
        # 			check_idx(fasta)
        # 		except:
        # 			print('No index for file: '+fasta+'. Please create!')
        # 		readnum = 0
        write_header(out)
        samfile = parse_bam(bam)
        header = ""

        try:
            samfile.check_index()
        except:
            print("No index for file: " + bam + ". Please create!")
        try:
            header = samfile.header["SQ"]

        except:
            print(
                "Could not read header of file "
                + bam
                + ". Multithreading only per file not per chromosome."
            )

        if "SQ" in header:
            chroms = []
            for k, v in header["SQ"].items():
                print(k, v)
                chroms.extend = str.split(":", k)[1]
            for chrom in chroms:
                pool.apply_async(
                    collect, args=(bam, fasta, out, chrom)
                )  # , callback = log_result)
        else:
            # 			res = pool.apply_async(test,args=(27,))
            res = pool.apply_async(
                collect,
                args=(
                    bam,
                    fasta,
                    out,
                ),
            )  # , callback = log_result)
    # 			collect(samfile, fasta, out)

    pool.close()
    pool.join()


# 	print(res.get())


def test(x):
    try:
        with open("bla", "a") as h:
            print(x ** 2, file=h)
    # 		raise Exception
    # 		return a
    except Exception as error:
        # capture the exception and bundle the traceback
        # in a string, then raise new exception with the string traceback
        raise Exception("".join(tb.format_exception(*sys.exc_info())))


def collect(bam, fasta, out, chrom=None):

    # 	with open('bla','a') as h:
    # 		print('Collected', file=h)

    statistics = collections.OrderedDict()
    statistics["reads"] = collections.OrderedDict()
    statistics["ref"] = collections.OrderedDict()

    samfile = parse_bam(bam)

    if chrom is not "None":
        for read in samfile.fetch(chrom):
            if not read.is_unmapped:
                chrom = read.reference_name
                cigar = read.cigarstring
                if "*" in cigar:
                    continue
                seq = read.query_alignment_sequence
                alignmentstart = read.reference_start

                if (len(statistics["reads"]) > 0) and (
                    chrom not in statistics["reads"]
                ):  # or (alignmentstart not in statistics['reads'][chrom])):
                    print("Fetching reference sequence and printing")

                    write_stats(out, statistics, fasta)
                    statistics["reads"].clear()
                    statistics["ref"].clear()

                if chrom not in statistics["reads"]:
                    print("Parsing Chromosome " + chrom)

                get_stats(chrom, cigar, seq, alignmentstart, statistics)

    else:
        for read in samfile.fetch():
            if not read.is_unmapped:
                chrom = read.reference_name
                cigar = read.cigarstring
                if "*" in cigar:
                    continue
                seq = read.query_alignment_sequence
                alignmentstart = read.reference_start

                if (len(statistics["reads"]) > 0) and (
                    chrom not in statistics["reads"]
                ):  # or (alignmentstart not in statistics['reads'][chrom])):
                    print("Fetching reference sequence and printing")

                    write_stats(out, statistics, fasta)
                    statistics["reads"].clear()
                    statistics["ref"].clear()

                if chrom not in statistics["reads"]:
                    print("Parsing Chromosome " + chrom)

                get_stats(chrom, cigar, seq, alignmentstart, statistics)

    write_stats(out, statistics, fasta)
    close_bam(samfile)


# 	print('Reads parsed: ' + str(readnum))


def parse_bam(bam):
    if ".bam" in bam:
        return pysam.AlignmentFile(bam, "rb")
    elif ".sam" in bam:
        return pysam.AlignmentFile(bam, "r")


def read_head(bam):
    if ".bam" in bam:
        return pysam.AlignmentFile(bam, "rb").header
    elif ".sam" in bam:
        return pysam.AlignmentFile(bam, "r").header


def write_bam(bam, template):
    return pysam.AlignmentFile(bam, "wh", template=template)


def close_bam(samfile):
    samfile.close()


def check_idx(file):
    if (
        (os.path.isfile(file + ".idx"))
        or (os.path.isfile(str.join(".", str.split(".", file)[0:-1], "fai")))
        or (os.path.isfile(str.join(".", str.split(".", file)[0:-1], "faidx")))
    ):
        return True
    else:
        return False


def get_stats(chrom, cigar, seq, alignmentstart, statistics):

    # 	with open('bla','a') as h:
    # 		print('Fetching Stats', file=h)

    pattern = re.compile("([0-9]*)([DMINSHP=XB])")
    pos = alignmentstart + 1
    char = 0
    if chrom not in statistics["reads"]:
        statistics["reads"][chrom] = collections.OrderedDict()

    for n, c in pattern.findall(cigar):
        if n:
            n = int(n)
        else:
            n = 1

        if c == "M" or c == "X" or c == "=":
            for i in range(n):
                if pos not in statistics["reads"][chrom]:
                    statistics["reads"][chrom][pos] = collections.OrderedDict()
                if seq[char] not in statistics["reads"][chrom][pos]:
                    statistics["reads"][chrom][pos][seq[char]] = 0

                statistics["reads"][chrom][pos][seq[char]] += 1
                char += 1
                pos += 1
        if c in ["D", "N"]:
            pos += n
        if c in ["I"]:
            char += n
        if c in ["H", "P", "B", "S"]:  # No idea what PB mean
            continue


# 	biggestseen = pos


def get_ref(fa, statistics):
    faseq = Fasta(fa)

    # 	with open('bla','a') as h:
    # 		print('Fetching Ref', file=h)

    for chrom in statistics["reads"]:
        if chrom not in statistics["ref"]:
            statistics["ref"][chrom] = collections.OrderedDict()
        for pos in statistics["reads"][chrom]:
            if pos not in statistics["ref"][chrom]:
                statistics["ref"][chrom][pos] = collections.OrderedDict()
            nuc = faseq[chrom][pos - 1].seq.upper()
            statistics["ref"][chrom][pos] = nuc


def write_header(bam):
    head = (
        str.join(
            "\t",
            ["Chromosome", "Position", "Ref", "A", "C", "G", "T", "N", "RelRefFreq"],
        )
        + "\n"
    )
    o = gzip.open(bam, "ab")
    o.write(bytes(head, encoding="UTF-8"))
    o.close()


def write_stats(bam, statistics, fasta, start=None):

    o = gzip.open(bam, "ab")
    outstr = ""
    coverage = ""

    get_ref(fasta, statistics)

    for chrom in statistics["reads"]:
        for pos in statistics["reads"][chrom]:
            if start is not None:
                if int(pos) > int(start):
                    continue
            outstr = (
                str.join("\t", [str(chrom), str(pos), statistics["ref"][chrom][pos]])
                + "\t"
            )
            sum = 0
            for nuc in ["A", "C", "G", "T", "N"]:
                if nuc in statistics["reads"][chrom][pos]:
                    sum += statistics["reads"][chrom][pos][nuc]
                    coverage += str(statistics["reads"][chrom][pos][nuc]) + "\t"
                else:
                    coverage += str(0) + "\t"
            if statistics["ref"][chrom][pos] in statistics["reads"][chrom][pos]:
                coverage += str(
                    statistics["reads"][chrom][pos][statistics["ref"][chrom][pos]] / sum
                )
            else:
                coverage += str(0)
            outstr += coverage + "\n"
            o.write(bytes(outstr, encoding="UTF-8"))
            outstr = ""
            coverage = ""
    ###Must not change size during iteration
    # 			del statistics['reads'][chrom][pos]
    # 			del statistics['ref'][chrom][pos]
    o.close()


if __name__ == "__main__":
    args = parseargs()
    collectstats(args.fasta, args.bams, args.outdir, args.procs, args.verbosity)

#
# CollectBamStat.py ends here
