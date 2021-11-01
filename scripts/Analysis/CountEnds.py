# CountEnds.py ---
#
# Filename: CollectBamStat.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue May 15 12:44:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Jan 30 15:53:25 2019 (+0100)
#           By: Joerg Fallmann
#     Update #: 923
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
import logging

# import importlib
from multiprocessing import Pool, Manager
import traceback as tb

# import pickle
# Container
import collections
import inspect

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


### MAIN
def parseargs():
    parser = argparse.ArgumentParser(
        description="Read BAM file read by read and collect mapping statistics"
    )
    parser.add_argument(
        "-f", "--fasta", type=str, default=None, help="Fasta with cluster info"
    )
    parser.add_argument(
        "-b",
        "--bams",
        type=str,
        help="Mapped reads BAM(s) comma separated, need to besorted and indexed",
    )
    parser.add_argument(
        "-l",
        "--lookup",
        type=str,
        default=None,
        help="Search for given ID in chrom tag of reads, e.g. Cluster, chr1, etc.",
    )
    parser.add_argument(
        "-t", "--offset", type=int, default=0, help="Offset for search around read end"
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


def collectends(bams, outdir, procs, verbosity, offset, fasta=None, lookup=None):
    try:
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
        pool = Pool(processes=num_processes, maxtasksperchild=1)
        collector = collections.OrderedDict()
        translater = None

        if fasta:
            translater = cluster2trna(fasta)[
                "cluster"
            ]  # Make sure to use the right fasta for the bams

        for bam in bams.split(","):
            fn = os.path.basename(bam)
            out = os.path.join(outdir, fn + ".ends.gz")
            collector[fn] = collections.OrderedDict()
            collector[fn]["out"] = out
            collector[fn]["res"] = collections.OrderedDict()

            if os.path.isfile(out):
                print("File " + out + " exists, will be deleted")
                os.remove(out)

            write_header(out, translater)

            if os.path.exists(bam) and os.path.getsize(bam) > 0:
                samfile = parse_bam(bam)
                try:
                    ####Check if index files available
                    samfile.check_index()
                except:
                    print("No index for file: " + bam + ". Please create!")
                    sys.exit()

                header = None
                try:
                    if "SQ" in samfile.header:
                        header = samfile.header["SQ"]
                except:
                    print(
                        "Could not read header of file "
                        + bam
                        + ". Multithreading only per file not per chromosome."
                    )

                if header:
                    for i, chrom in enumerate(h["SN"] for h in header):
                        if lookup:
                            if lookup not in chrom:
                                continue
                        res = pool.apply_async(
                            collect, args=(bam, offset, chrom), error_callback=eprint
                        )
                        rescheck = res.get()
                        if chrom in rescheck:
                            collector[fn]["res"][chrom] = rescheck[chrom]

                else:
                    res = pool.apply_async(
                        collect, args=(bam, offset), error_callback=eprint
                    )
                    rescheck = res.get()
                    if len(rescheck) > 0:
                        collector[fn]["res"] = rescheck
            else:
                write_empty(out)
        pool.close()
        pool.join()

        for bam in collector:
            write_stats(collector[bam], translater)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def collect(bam, offset, chrom=None):
    # READS HAVE TO MAP TO END OF CLUSTER
    # TODO: CHECK IF NUMBERS ARE CORRECT!
    try:
        stats = collections.OrderedDict()

        samfile = parse_bam(bam)
        refs = dict(zip(samfile.references, samfile.lengths))

        if chrom and chrom is not "None":
            print("Parsing Chromosome parallel " + str(chrom))
            for read in samfile.fetch(chrom):
                #                if not read.is_unmapped:                        # This currently doesn't work with segemehl!
                seq = read.query_alignment_sequence.upper()  # read without softclip
                rawseq = read.query_sequence.upper()  # read with softclip
                refseq = (
                    read.get_reference_sequence().upper()
                )  # reference part of alignment

                if read.reference_end in range(
                    refs[chrom] - offset, refs[chrom] - offset + 6 + 1
                ):  # Should make sure we are at the read end but also allow softclips
                    tags = dict(read.get_tags())
                    alignmentend = str(seq[-6:])
                    rawend = str(rawseq[-6:])
                    refend = str(refseq[-6:])

                    if not chrom in stats:
                        stats[chrom] = collections.OrderedDict()
                        stats[chrom]["seq"] = collections.OrderedDict()
                        stats[chrom]["raw"] = collections.OrderedDict()
                        stats[chrom]["ref"] = collections.OrderedDict()

                    if (
                        alignmentend in stats[chrom]["seq"]
                    ):  # Collect ends and count them
                        stats[chrom]["seq"][alignmentend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["seq"][alignmentend] = 1 / tags["NH"]
                    if rawend in stats[chrom]["raw"]:  # Collect ends and count them
                        stats[chrom]["raw"][rawend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["raw"][rawend] = 1 / tags["NH"]
                    if refend in stats[chrom]["ref"]:  # Collect ends and count them
                        stats[chrom]["ref"][refend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["ref"][refend] = 1 / tags["NH"]

        else:
            for read in samfile.fetch():
                #                if not read.is_unmapped:  # This currently doesn's work with segemehl!
                chrom = read.reference_name
                seq = read.query_alignment_sequence.upper()
                rawseq = read.query_sequence.upper()
                refseq = (
                    read.get_reference_sequence().upper()
                )  # reference part of alignment

                if read.reference_end in range(
                    refs[chrom] - offset, refs[chrom] - offset + 6 + 1
                ):  # Should make sure we are at the read end but also allow softclips
                    #                    if refs[chrom]-offset == read.reference_end:
                    tags = dict(read.get_tags())
                    alignmentend = str(seq[-6:])
                    rawend = str(rawseq[-6:])
                    refend = str(refseq[-6:])

                    if not chrom in stats:
                        stats[chrom] = collections.OrderedDict()
                        stats[chrom]["seq"] = collections.OrderedDict()
                        stats[chrom]["raw"] = collections.OrderedDict()
                        stats[chrom]["ref"] = collections.OrderedDict()

                    if (
                        alignmentend in stats[chrom]["seq"]
                    ):  # Collect ends and count them
                        stats[chrom]["seq"][alignmentend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["seq"][alignmentend] = 1 / tags["NH"]
                    if rawend in stats[chrom]["raw"]:  # Collect ends and count them
                        stats[chrom]["raw"][rawend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["raw"][rawend] = 1 / tags["NH"]
                    if refend in stats[chrom]["ref"]:  # Collect ends and count them
                        stats[chrom]["ref"][refend] += 1 / tags["NH"]
                    else:
                        stats[chrom]["ref"][refend] = 1 / tags["NH"]

        close_bam(samfile)

        return stats

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def fillre(col):
    try:
        rc = col.get()
        return rc

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def collect_collector(bam, chrom=None):

    try:
        rc = readcollector(os.path.basename(bam))
        samfile = parse_bam(bam)
        if chrom and chrom is not "None":
            printlog("Parsing Chromosome parallel" + str(chrom))
            for read in samfile.fetch(chrom):
                if not read.is_unmapped:
                    seq = read.query_alignment_sequence
                    tags = dict(read.get_tags())
                    alignmentend = str(seq[-6:])
                    if not rc.checkchrom(chrom):
                        rc.add(chrom)
                    if rc.checkend(
                        [chrom, alignmentend]
                    ):  # Collect ends and count them
                        rc.addval([chrom, alignmentend], 1 / tags["NH"])
                    else:
                        rc.add([chrom, alignmentend], 1 / tags["NH"])
        #                        rc.addend(chrom, alignmentend, 1/tags['NH'])

        else:
            for read in samfile.fetch():
                if not read.is_unmapped:
                    chrom = read.reference_name
                    seq = read.query_alignment_sequence
                    tags = dict(read.get_tags())
                    alignmentend = seq[-6:]
                    if not rc.checkchrom(chrom):
                        rc.add(chrom)
                    if rc.checkend(
                        [chrom, alignmentend]
                    ):  # Collect ends and count them
                        rc.addval([chrom, alignmentend], 1 / tags["NH"])
                    else:
                        rc.add([chrom, alignmentend], 1 / tags["NH"])
        #           close_bam(samfile)
        return rc

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def parse_bam(bam):
    try:
        if ".bam" in bam:
            return pysam.AlignmentFile(bam, "rb")
        elif ".sam" in bam:
            if ".gz" in bam[-4:]:
                return pysam.AlignmentFile(gzip.open(bam, "rt"), "r")
            else:
                return pysam.AlignmentFile(bam, "r")
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def parse_cigar(cigarstring, length, start):

    pattern = re.compile("([0-9]*)([DMINSHP=XB])")
    stop = start

    for n, c in pattern.findall(cigar):
        if n:
            n = int(n)
        else:
            n = 1

        if c == "M" or c == "X" or c == "=":
            stop += n
        if c is "D":
            stop += n
            length += n
        if c is "N":
            stop += n

        if c in ["I"]:
            char += n
        if c in ["H", "P", "B", "S"]:  # No idea what PB mean
            continue

    return


def read_head(bam):
    try:
        if ".bam" in bam:
            return pysam.AlignmentFile(bam, "rb").header
        elif ".sam" in bam:
            return pysam.AlignmentFile(bam, "r").header
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_bam(bam, template):
    try:
        return pysam.AlignmentFile(bam, "wh", template=template)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def close_bam(samfile):
    samfile.close()


def check_idx(file):
    try:
        if (
            (os.path.isfile(file + ".idx"))
            or (os.path.isfile(str.join(".", split(".", file)[0:-1], "fai")))
            or (os.path.isfile(str.join(".", split(".", file)[0:-1], "faidx")))
        ):
            return True
        else:
            return False
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def get_stats(chrom, alignmentend, statistics):
    try:
        #   with open('bla','a') as h:
        #       print('Fetching Stats', file=h)

        if chrom not in statistics["ends"]:
            statistics["ends"][chrom] = collections.OrderedDict()

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

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def get_ref(fa, statistics):
    try:
        faseq = Fasta(fa)

        for chrom in statistics["reads"]:
            if chrom not in statistics["ref"]:
                statistics["ref"][chrom] = collections.OrderedDict()
            for pos in statistics["reads"][chrom]:
                if pos not in statistics["ref"][chrom]:
                    statistics["ref"][chrom][pos] = collections.OrderedDict()
                nuc = faseq[chrom][pos - 1].seq.upper()
                statistics["ref"][chrom][pos] = nuc
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_header(out, fastadict=None):
    try:
        if fastadict:
            head = str.join(
                "\t",
                [
                    "Chromosome",
                    "Seq",
                    "End",
                    "Count",
                    "ChromRelRefFreq",
                    "TotalRelRefFreq",
                    "tRNAs",
                ],
            )
        else:
            head = str.join(
                "\t",
                [
                    "Chromosome",
                    "Seq",
                    "End",
                    "Count",
                    "ChromRelRefFreq",
                    "TotalRelRefFreq",
                ],
            )

        with gzip.open(out, "wb") as o:
            o.write(bytes(head + "\n", encoding="UTF-8"))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_empty(out):
    try:
        with gzip.open(out, "ab") as o:
            o.write(bytes("0\n", encoding="UTF-8"))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def write_stats(stat, fastadict=None):
    try:
        coverage = ""
        allends = collections.OrderedDict()
        totalcount = {}
        readnum = {}
        if stat is not None:
            out = stat.pop("out", None)  # get key and remove from dict
            statistics = stat.pop("res", None)  # get key and remove from dict

            for chrom in statistics:
                if chrom not in allends:
                    allends[chrom] = collections.OrderedDict()
                    allends[chrom]["seq"] = collections.OrderedDict()
                    allends[chrom]["raw"] = collections.OrderedDict()
                    allends[chrom]["ref"] = collections.OrderedDict()
                if chrom not in totalcount:
                    totalcount[chrom] = collections.OrderedDict()
                    totalcount[chrom]["seq"] = 0
                    totalcount[chrom]["raw"] = 0
                    totalcount[chrom]["ref"] = 0

                for rtype in statistics[chrom]:
                    if not rtype in readnum:
                        readnum[rtype] = 0
                    for end in statistics[chrom][rtype]:
                        totalcount[chrom][rtype] += statistics[chrom][rtype][end]
                        readnum[rtype] += statistics[chrom][rtype][end]
                        allends[chrom][rtype][end] = statistics[chrom][rtype][end]
                        for x in range(-5, 0, 1):
                            newend = end[x:]
                            if newend in allends[chrom][rtype]:
                                allends[chrom][rtype][newend] += statistics[chrom][
                                    rtype
                                ][end]
                            else:
                                allends[chrom][rtype][newend] = statistics[chrom][
                                    rtype
                                ][end]

            outstr = []
            for chrom in allends:
                for rtype in allends[chrom]:
                    for end in allends[chrom][rtype]:
                        if allends[chrom][rtype][end] is not "":
                            tosave = None
                            if fastadict:
                                info = ""
                                for id in fastadict:
                                    if (id.upper() == chrom.upper()) or any(
                                        [
                                            x.upper() in chrom.upper().split(".")[-1]
                                            for x in fastadict[id]
                                        ]
                                    ):
                                        info = fastadict[id]
                                tosave = str.join(
                                    "\t",
                                    [
                                        str(chrom),
                                        rtype,
                                        end,
                                        str(allends[chrom][rtype][end]),
                                        str(
                                            allends[chrom][rtype][end]
                                            / totalcount[chrom][rtype]
                                        ),
                                        str(
                                            allends[chrom][rtype][end] / readnum[rtype]
                                        ),
                                        str(",".join(info)),
                                    ],
                                )
                            else:
                                tosave = str.join(
                                    "\t",
                                    [
                                        str(chrom),
                                        rtype,
                                        end,
                                        str(allends[chrom][rtype][end]),
                                        str(
                                            allends[chrom][rtype][end]
                                            / totalcount[chrom][rtype]
                                        ),
                                        str(
                                            allends[chrom][rtype][end] / readnum[rtype]
                                        ),
                                    ],
                                )
                            if tosave and len(tosave) > 1:
                                outstr.append(tosave)

            with gzip.open(out, "ab") as o:
                o.write(bytes("\n".join(outstr), encoding="UTF-8"))

        else:
            with gzip.open(out, "wb") as o:
                o.write(bytes("", encoding="UTF-8"))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


def printlog(msg):
    with open("log", "a") as l:
        print(str(msg), file=l)


if __name__ == "__main__":
    args = parseargs()
    collectends(
        args.bams,
        args.outdir,
        args.procs,
        args.verbosity,
        args.offset,
        args.fasta,
        args.lookup,
    )

#
# CountEnds.py ends here
