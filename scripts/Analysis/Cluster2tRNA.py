# Cluster2tRNA.py ---
#
# Filename: Cluster2tRNA.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Oct  9 15:18:52 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Oct  9 21:36:52 2018 (+0200)
#           By: Joerg Fallmann
#     Update #: 40
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
from io import StringIO
import os
import gzip
import re
import sys
import string
import collections
import traceback as tb
import inspect

####load own modules
cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../../MONSDA",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
# sys.path=[str(os.getenv('HOME')+"/Work/Scripts/Python/lib")] + sys.path
from Collection import *

####Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq

####numpy and matplolib and pyplot
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
from Logger import *

### MAIN
def parseargs():
    parser = argparse.ArgumentParser(
        description="Read BAM file read by read and collect mapping statistics"
    )
    parser.add_argument(
        "-f", "--fasta", type=str, help="Reference genome FASTA, needs to be indexed"
    )
    parser.add_argument(
        "-o", "--outdir", type=str, default=".", help="Output directory"
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


if __name__ == "__main__":
    try:
        args = parseargs()
        if args.outdir:
            print("Checking or creating outdir " + str(args.outdir))
            if not os.path.isabs(args.outdir):
                outdir = os.path.abspath(args.outdir)
            if not os.path.exists(args.outdir):
                os.makedirs(args.outdir)
        else:
            outdir = os.path.abspath(os.getcwd())

        translater = cluster2trna(args.fasta)
        outstring = []

        for cluster in translater["cluster"]:

            outstring.append(
                "\t".join([cluster, ",".join(translater["cluster"][cluster])])
            )

        for chrom in translater["tRNA"]:
            for strand in translater["tRNA"][chrom]:
                for tRNA in translater["tRNA"][chrom][strand]:
                    outstring.append("\t".join([tRNA, chrom, strand]))

        with open(os.path.join(outdir, "cluster2trna_" + args.fasta), "w") as o:
            print("\n".join(outstring), file=o)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type,
            exc_value,
            exc_tb,
        )
        with open("error", "a") as h:
            print("".join(tbe.format()), file=h)


#
# Cluster2tRNA.py ends here
