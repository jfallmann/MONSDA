#!/usr/bin/env python3

import sys
import argparse
import os
import inspect
import gzip
import glob
import re
import logging
import traceback as tb

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
        description="Parses FeatureCounts output tables into DIEGO ready format"
    )
    # parser.add_argument("-l", "--list", type=str, required=True, help="List of samples")
    parser.add_argument(
        "-n",
        "--sample_name",
        action="store_true",
        help=" provide -n if sample names instead of group names should be used for header",
    )
    parser.add_argument(
        "-o",
        "--order",
        action="store_true",
        help="if wanted the order of conditions can be given as comma separated list",
    )
    parser.add_argument(
        "-c", "--conditions", required=True, type=str, help="Conditions to compare"
    )
    parser.add_argument(
        "-t", "--types", required=False, type=str, help="Sequencing types to compare"
    )
    parser.add_argument(
        "-b", "--batches", required=False, type=str, help="Sample batches to compare"
    )
    parser.add_argument(
        "-r",
        "--replicates",
        required=True,
        type=str,
        help="Replicates belonging to conditions",
    )
    parser.add_argument(
        "--cutoff", dest="cutoff", type=int, default=0, help="cutoff for minimum count"
    )
    parser.add_argument(
        "-p",
        "--paired",
        required=False,
        type=str,
        default=None,
        help="Sequencing strategy for sample name processing",
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
    parser.add_argument("--loglevel", default="DEBUG", help="Log verbosity")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


class Sample_list(object):

    group_name = ""
    replicate_names = []
    replicate_paths = []
    # the class constructor
    def __init__(self, group_name):
        self.group_name = group_name
        self.replicate_names = list()
        self.replicate_paths = list()
        self.replicate_types = list()
        self.replicate_batches = list()


def prepare_table(
    conditions,
    replicates,
    types,
    batches,
    paired,
    table,
    anno,
    sample_name=None,
    order=None,
    cutoff=None,
):
    try:
        # slist,
        logid = scriptname + ".prepare_table: "
        my_groups = {}
        list_size = 0

        # CLEANUP
        oldtab = os.path.abspath(table)
        oldanno = os.path.abspath(anno)
        for oldfile in glob.glob(oldtab):
            os.rename(oldfile, oldfile + ".bak")
            log.warning(
                logid
                + "Found old table file"
                + oldfile
                + ", was moved to "
                + oldfile
                + ".bak"
            )
        for oldfile in glob.glob(oldanno):
            os.rename(oldfile, oldfile + ".bak")
            log.warning(
                logid
                + "Found old anno file"
                + oldfile
                + ", was moved to "
                + oldfile
                + ".bak"
            )

        replist = str(replicates).strip().split(",")
        condlist = str(conditions).strip().split(",")
        typelist = str(types).strip().split(",") if types is not None else None
        batchlist = str(batches).strip().split(",") if batches is not None else None
        pairedlist = str(paired).strip().split(",") if paired is not None else None

        log.debug(logid + "REPS: " + str(replist) + "\tLEN: " + str(len(replist)))
        log.debug(logid + "CONDS: " + str(condlist) + "\tLEN: " + str(len(condlist)))

        if types is not None:
            log.debug(
                logid + "TYPES: " + str(typelist) + "\tLEN: " + str(len(typelist))
            )

        if batches is not None:
            log.debug(
                logid + "BATCHES: " + str(batchlist) + "\tLEN: " + str(len(batchlist))
            )

        if paired is not None:
            log.debug(
                logid + "PAIRED: " + str(pairedlist) + "\tLEN: " + str(len(pairedlist))
            )

        for i in range(len(replist)):
            rep = None
            cond = None
            typ = None
            bat = None

            rep = str(replist[i])
            cond = str(condlist[i])
            typ = str(typelist[i]) if types is not None else None
            bat = str(batchlist[i]) if batches is not None else None

            if not rep or not cond:
                log.warning(logid + "No rep/cond found for sample " + str(replist[i]))

            log.debug(logid + "rep/cond: " + str([rep, cond]))

            list_size += 1

            if cond in my_groups:
                my_groups[cond].replicate_paths.append(rep)
                my_groups[cond].replicate_names.append(str.split(os.sep, rep)[-1])
                if typ is not None:
                    my_groups[cond].replicate_types.append(typ)
                if bat is not None:
                    my_groups[cond].replicate_batches.append(bat)
            else:
                my_groups[cond] = make_sample_list(cond)
                my_groups[cond].replicate_paths.append(rep)
                my_groups[cond].replicate_names.append(str.split(os.sep, rep)[-1])
                if typ is not None:
                    my_groups[cond].replicate_types.append(typ)
                if bat is not None:
                    my_groups[cond].replicate_batches.append(bat)

        log.debug(logid + "MyGroups: " + str(my_groups.keys()))

        myMatrix = []
        myMatrix.append([])
        myMatrix[0].append("\t".join(["junction", "type"]))
        sample_counter = 0

        conds = []
        if order:
            conds = order.split(",")
        else:
            conds = [x for x in my_groups.keys()]

        log.debug(logid + "CONDS: " + str(conds))
        typeanno = list()
        batchanno = list()
        for gruppies in conds:
            condition_index = -1
            rep_nr = 0
            for replicates in my_groups[gruppies].replicate_paths:
                log.info(logid + "Processing: " + str(replicates))
                condition_index += 1
                sample_counter += 1
                rep_nr += 1

                if sample_name:
                    myMatrix[0].append(
                        my_groups[gruppies].replicate_names[condition_index]
                    )
                    typeanno.append(
                        my_groups[gruppies].replicate_types[condition_index]
                    )
                    batchanno.append(
                        my_groups[gruppies].replicate_batches[condition_index]
                    )
                else:
                    myMatrix[0].append(
                        str(my_groups[gruppies].group_name) + "_" + str(rep_nr)
                    )
                    typeanno.append(
                        my_groups[gruppies].replicate_types[condition_index]
                    )
                    batchanno.append(
                        my_groups[gruppies].replicate_batches[condition_index]
                    )
                if ".gz" in replicates:
                    myInput = gzip.open(replicates, "rt")
                else:
                    myInput = open(replicates, "r")

                lineNumber = 0
                oldgene = ""
                junctionnr = 0
                for line in myInput:
                    if "#" in line[0:5] or ".bam" in line[-10:]:
                        continue
                    columns = line.strip().split("\t")
                    if columns[0] != "name" and columns[1] != "count":
                        lineNumber += 1
                        gene = str(columns[0]).replace(":", "|")
                        if gene != oldgene:
                            junctionnr = 0
                            oldgene = gene
                        else:
                            junctionnr += 1
                        chr = (
                            str(columns[1]).replace(":", "|")
                            + "-"
                            + str(columns[2])
                            + "-"
                            + str(columns[3])
                        )
                        id = (
                            str(columns[4]).replace(":", "|")
                            + "junction_"
                            + str(junctionnr)
                        )
                        if sample_counter == 1:
                            newListi = []
                            myMatrix.append(newListi)
                            myMatrix[lineNumber].append(":".join([gene, id, chr]))
                        if len(columns) > 1 and columns[-1] != columns[0]:
                            myMatrix[lineNumber].append(round(float(columns[-1])))
                        else:
                            myMatrix[lineNumber].append("0")

        line = "\t".join(["\t".join(myMatrix[0]), "geneID", "geneName"])
        annos = list()

        for i in range(1, len(myMatrix[0])):
            c = myMatrix[0][i]
            a = str.join("_", str(c).split("_")[:-1])
            a += "\t" + str(typeanno[i - 1]) if types is not None else None
            a += "\t" + str(batchanno[i - 1]) if batches is not None else None
            annos.append(str(c) + "\t" + str(a))

        with gzip.open(anno, "wb") as a:
            a.write(bytes("\n".join(annos) + "\n", encoding="UTF-8"))

        toprint = str(line) + "\n"
        for z in range(1, len(myMatrix)):
            zeilen = myMatrix[z]
            if "geneName" in str(zeilen):
                continue
            willprint = False
            gene, id, chr = str(zeilen[0]).split(":")
            line = "\t".join([str(chr).replace("-", ":", 1) + ":" + id, "N_w" + "\t"])
            for x in range(1, len(zeilen)):
                line = line + str(zeilen[x]) + "\t"
                if int(zeilen[x]) >= cutoff:
                    willprint = True
            line += "\t".join([gene, gene])
            if willprint:
                toprint = toprint + str(line) + "\n"

        with gzip.open(table, "wb") as t:
            t.write(bytes(str(toprint), encoding="UTF8"))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))


def make_sample_list(group_name):
    sample_list = Sample_list(group_name)
    return sample_list


####################
####    MAIN    ####
####################

if __name__ == "__main__":

    logid = scriptname + ".main: "
    try:
        args = parseargs()
        try:
            makelogdir("LOGS")
            log = setup_logger(
                name=scriptname,
                log_file="LOGS/" + scriptname + ".log",
                logformat="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
                datefmt="%m-%d %H:%M",
                level=args.loglevel,
            )
            log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        prepare_table(
            args.conditions,
            args.replicates,
            args.types,
            args.batches,
            args.paired,
            args.table,
            args.anno,
            args.sample_name,
            args.order,
            args.cutoff,
        )
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        log.error(logid + "".join(tbe.format()))

#
# FeatureCounts2DIEGO.py ends here
