# GenerateTrackDb.py ---
#
# Filename: GenerateTrackDb.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Dec  4 09:54:46 2017 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Mar  9 10:23:08 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 161
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
###############
###Imports
###############
import sys
import os
import glob
import time
import gzip
import string
import argparse
import random
import traceback as tb
import inspect
import logging

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

###############
###Variables
###############

verbose = 0
date = time.strftime("%d/%m/%Y")

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

###############
###Command Line Options
###############


def parseargs():
    parser = argparse.ArgumentParser(description="GenerateTrackDb.py")
    parser.add_argument(
        "-e",
        "--header",
        type=bool,
        default=False,
        help="If header should be created or just entries added to the trackdb file",
    )
    parser.add_argument(
        "-t",
        "--track",
        type=str,
        default="trackdb.txt",
        help="Name of the trackdb file to generate or to extend",
    )
    parser.add_argument(
        "-f",
        "--filenames",
        type=str,
        help="Name of the file containing bwfilenames that should be added to the trackdb file, can be STDIN",
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="Autogentrack",
        help="Name of the track to create, can be anything",
    )
    parser.add_argument(
        "-s",
        "--shortlabel",
        type=str,
        help="Short label for track and tracks, default is first 20chars of bw files",
    )
    parser.add_argument(
        "-l",
        "--longlabel",
        type=str,
        help="Long label for track and tracks, default is name of bw files",
    )
    parser.add_argument(
        "-u",
        "--baseurl",
        type=str,
        help="Base url that will be used for all the tracks",
    )
    parser.add_argument(
        "-g",
        "--genome",
        type=str,
        help="Genome used for analysis, will be added to genomes.txt",
    )
    parser.add_argument(
        "-b",
        "--hub",
        type=str,
        default="AutogenHub",
        help="Will generate hub.txt file.",
    )
    parser.add_argument(
        "-m", "--mail", type=str, default="is@egal.com", help="Email of track owner."
    )
    parser.add_argument(
        "-i",
        "--uid",
        type=str,
        default="",
        help="Unique ID for output in case of bulk generation",
    )
    parser.add_argument(
        "-x",
        "--splitbw",
        type=bool,
        default=False,
        help="Will split into reverse and forward tracks, default is to negate reverse strand values.",
    )
    return parser.parse_args()


###############
###MAIN
###############


def main(
    header,
    track,
    shortlabel,
    longlabel,
    filenames,
    name,
    url,
    genome,
    hub,
    mail,
    splitbw,
    uid,
):
    try:
        if header and splitbw:
            if not shortlabel:
                shortlabel = os.path.split(name)[1].split(".")[0]
            if not longlabel:
                longlabel = name

            head = "\n".join(
                [
                    "# TrackDb_generated_by_GenerateTrackDb.txt",
                    "track {tr}",
                    "container multiWig",
                    "noInherit on",
                    "shortLabel {sl}",
                    "longLabel {ll}",
                    "type bigWig",
                    "configurable on",
                    "visibility full",
                    "aggregate transparentOverlay",
                    "showSubtrackColorOnUi on",
                    "autoScale on",
                    "windowingFunction maximum",
                    "priority 1000",
                    "alwaysZero on",
                    "yLineMark 0",
                    "yLineOnOff on",
                    "maxHeightPixels 200:64:32",
                ]
            ).format(tr=name, sl=shortlabel, ll=longlabel)

            with open(uid + "re_" + track, "a") as out:
                out.write(head + "\n\n")
            with open(uid + "fw_" + track, "a") as out:
                out.write(head + "\n\n")

        if filenames == "STDIN" or filenames == "-":
            bigwigs = []
            for line in sys.stdin:
                bigwigs.append(line.split("\n", 2)[0])

        elif os.path.isfile(filenames):
            bigwigs = []
            if ".gz" in filenames:
                b = gzip.open(filenames, "rb")
                for line in b:
                    bigwigs.append(line.split("\n", 2)[0])
                    b.close()
            else:
                b = open(filenames, "r")
                for line in b:
                    bigwigs.append(line.split("\n", 2)[0])
                    b.close()

        fwtracks = []
        retracks = []
        head = ""
        idx = 1

        for content in bigwigs:
            bw = content
            sample = os.path.split(bw)[1].split(".")[0]
            colorfw = "128,0,0"  # random_color()
            colorre = "0,128,0"
            shortlabel = sample
            longlabel = bw.split(".bw", 2)[0]

            if "\t" in bw:
                bw, sample = content.split("\t", 2)

            if header:
                if not splitbw and ".re" in longlabel:
                    content = "\n".join(
                        [
                            "track " + sample + "_re",
                            "bigDataUrl {1}/" + bw,
                            "shortLabel {2}",
                            "longLabel {3}",
                            "type bigWig",
                            "negateValues on #uncomment of this is wanted",
                            "parent {4}",
                            "color {5}",
                        ]
                    ).format(bw, url, shortlabel, longlabel, sample, colorre)
                elif not splitbw and not ".re" in longlabel:
                    head = "\n".join(
                        [
                            "# TrackDb_generated_by_GenerateTrackDb.txt",
                            "track {tr}",
                            "container multiWig",
                            "noInherit on",
                            "shortLabel {sl}",
                            "longLabel {ll}",
                            "type bigWig",
                            "configurable on",
                            "visibility full",
                            "aggregate transparentOverlay",
                            "showSubtrackColorOnUi on",
                            "autoScale on",
                            "windowingFunction maximum",
                            "priority 1000",
                            "alwaysZero on",
                            "yLineMark 0",
                            "yLineOnOff on",
                            "maxHeightPixels 200:64:32",
                        ]
                    ).format(tr=sample, sl=shortlabel, ll=longlabel)

                    content = "\n".join(
                        [
                            "track " + sample + "_fw",
                            "bigDataUrl {1}/" + bw,
                            "shortLabel {2}",
                            "longLabel {3}",
                            "type bigWig",
                            "#negateValues on #uncomment of this is wanted",
                            "parent {4}",
                            "color {5}",
                        ]
                    ).format(bw, url, shortlabel, longlabel, sample, colorfw)
            else:
                if not splitbw and ".re" in longlabel:
                    content = "\n".join(
                        [
                            "track " + sample + "_re" + str(idx),
                            "bigDataUrl {1}/" + bw,
                            "shortLabel {2}" + str(idx),
                            "longLabel {3}" + str(idx),
                            "type bigWig",
                            "negateValues on #uncomment of this is wanted",
                            "color {5}",
                        ]
                    ).format(bw, url, shortlabel, longlabel, name, colorre)
                else:
                    content = "\n".join(
                        [
                            "track " + sample + "_fw" + str(idx),
                            "bigDataUrl {1}/" + bw,
                            "shortLabel {2}" + str(idx),
                            "longLabel {3}" + str(idx),
                            "type bigWig",
                            "#negateValues on #uncomment of this is wanted",
                            "color {5}",
                        ]
                    ).format(bw, url, shortlabel, longlabel, name, colorfw)

            idx += 1

            if splitbw and ".re" in longlabel:
                retracks.append(content)
            elif ".re" in longlabel:
                fwtracks.append(content)
            else:
                fwtracks.append(head + "\n\n" + content)

        if len(retracks) > 0:
            for i in retracks:
                with open(uid + "re_" + track, "a") as out:
                    out.write(i + "\n\n")

            for i in fwtracks:
                with open(uid + "fw_" + track, "a") as out:
                    out.write(i + "\n\n")

        else:
            for i in fwtracks:
                with open(uid + track, "a") as out:
                    out.write(i + "\n\n")

        if genome:
            if retracks:
                gen = "\n".join(
                    ["genome " + genome, "trackDb " + uid + "fw_{0}"]
                ).format(track)
                with open(uid + "genomes_fw.txt", "w") as out:
                    out.write(gen + "\n")

                gen = "\n".join(
                    ["genome " + genome, "trackDb " + uid + "re_{0}"]
                ).format(track)
                with open(uid + "genomes_re.txt", "w") as out:
                    out.write(gen + "\n")

            else:
                gen = "\n".join(["genome " + genome, "trackDb " + uid + "{0}"]).format(
                    track
                )
                with open(uid + "genomes.txt", "w") as out:
                    out.write(gen + "\n")

        if hub:
            if retracks:
                hubtxt = "\n".join(
                    [
                        "hub fw_" + hub,
                        "shortLabel FW_{0}",
                        "longLabel FW_{1}",
                        "genomesFile " + uid + "genomes_fw.txt",
                        "email {2}",
                    ]
                ).format(shortlabel, longlabel, mail)
                with open(uid + "hub_fw.txt", "w") as out:
                    out.write(hubtxt + "\n")

                hubtxt = "\n".join(
                    [
                        "hub re_" + hub,
                        "shortLabel RE_{0}",
                        "longLabel RE_{1}",
                        "genomesFile " + uid + "genomes_re.txt",
                        "email {2}",
                    ]
                ).format(shortlabel, longlabel, mail)
                with open(uid + "hub_re.txt", "w") as out:
                    out.write(hubtxt + "\n")

            else:
                hubtxt = "\n".join(
                    [
                        "hub " + hub,
                        "shortLabel {0}",
                        "longLabel {1}",
                        "genomesFile " + uid + "genomes.txt",
                        "email {2}",
                    ]
                ).format(shortlabel, longlabel, mail)
                with open(uid + "hub.txt", "w") as out:
                    out.write(hubtxt + "\n")
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        print("".join(tbe.format()))


###############
###Subs
###############
def random_color():
    rgbl = [255, 0, 0]
    random.shuffle(rgbl)
    return str(rgbl)


###############
###Standalone
###############

if __name__ == "__main__":
    args = parseargs()
    main(
        args.header,
        args.track,
        args.shortlabel,
        args.longlabel,
        args.filenames,
        args.name,
        args.baseurl,
        args.genome,
        args.hub,
        args.mail,
        args.splitbw,
        args.uid,
    )

##################################END################################


#
# GenerateTrackDb.py ends here
