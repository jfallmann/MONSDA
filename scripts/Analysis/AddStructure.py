#/usr/bin/env python3
# AddStructure.py ---
#
# Filename: AddStructure.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 10 18:00:42 2019 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Wed Sep 11 09:10:01 2019 (+0200)
#           By: Joerg Fallmann
#     Update #: 33
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
###Imports
import sys,os
import argparse
import traceback as tb
import gzip
import RNA

###Arguments
def parseargs():
	parser = argparse.ArgumentParser(description='Add structure to (bed) file containing sequence')
	parser.add_argument("-f", "--field", type=int, default=0, help='Which field contains the sequence, default is last field (0)')
	parser.add_argument("-b", "--bed", type=str, help='Bed or other tab separated file containing the sequence')

	return parser.parse_args()

###CODE

def addseq(field, bed):
    try:
        entries = parse_bed(bed)
        for line in entries:
            sequence = line.rstrip().split('\t')[field-1].upper()
            # create model details
            md = RNA.md()
            # create new fold_compound object
            fc = RNA.fold_compound(sequence, md)
            # compute minimum free energy (mfe) and corresponding structure
            (ss, mfe) = fc.mfe()

            sys.stdout.write('\t'.join([line.strip(),ss])+'\n')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def parse_bed(bed, annotated=None):
    try:
        if os.path.isfile(os.path.abspath(bed)):
            if '.gz' in bed:
                return gzip.open(bed,'rt')
            else:
                return open(bed,'rt')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


###MAIN
if __name__ == '__main__':
	args=parseargs()
	addseq(args.field, args.bed)

#
# AddStructure.py ends here
