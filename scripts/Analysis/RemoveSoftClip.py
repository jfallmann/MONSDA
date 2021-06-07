# RemoveSoftClip.py ---
#
# Filename: RemoveSoftClip.py
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
import pprint # pp = pprint.PrettyPrinter(indent=4) #pp.pprint(stuff)
from io import StringIO
import os
import gzip
import re
import sys
import pysam
from pyfaidx import Fasta
import traceback
#Container
import collections

### MAIN

def parseargs():
    parser = argparse.ArgumentParser(description='Read BAM file read by read and remove softclips')
    parser.add_argument("-f", "--fasta", type=str, help='Reference genome FASTA, needs to be indexed')
    parser.add_argument("-b", "--bams", type=str, help='Mapped reads BAM(s) comma separated, need to besorted and indexed')
    parser.add_argument("-o", "--outdir", type=str, default='.', help='Output directory')
    parser.add_argument("-c", "--cluster", type=str, default=None, help="Parses cluster information from chromosome tag")

    return parser.parse_args()

def process(fasta, bams, outdir, cluster=None):

    if outdir:
        print('Checking or creating outdir ' + str(outdir))
    if not os.path.isabs(outdir):
        outdir =  os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        outdir = os.path.abspath(os.getcwd())

    for bam in bams.split(','):
        fn = os.path.basename(bam)
        out = os.path.join(outdir, fn+'_nosoftclip.bam')

    if (os.path.isfile(out)):
        print('File '+ out +' exists, will be deleted')
        os.remove(out)

    samfile = parse_bam(bam)

    try:
        samfile.check_index()
    except:
        print('No index for file: '+bam+'. Please create!')

    remove_clip(bam, fasta, out, cluster)


def remove_clip(bam, fasta, out, cluster=None):

    samfile = parse_bam(bam)
    #header = samfile.header.as_dict()

    #if cluster:  # SQ tags need to be edited
    #    for k, v in header['SQ'].items():
    #        s, t, n, chrom, coord = k.split(':')
    #        k = chrom

    with pysam.Samfile(out, "wb", template=samfile) as outfile:

        for read in samfile.fetch():
            newread = pysam.AlignedSegment()
            if not read.is_unmapped:
                chrom = read.reference_name
                start, end = (0, 0)

                if cluster:  # Hammerhead_1::SM_V7_1:2251747-2251831(+)
                    t, n, chrom, coord = chrom.split(':')
                    start, end = map(int, coord.split('(')[0].split('-'))

                cigar = read.cigarstring
                if '*' in cigar:
                    continue

                newread.query_qualities = read.query_alignment_qualities
                newread.query_sequence = read.query_alignment_sequence
                newread.query_length = read.query_alignment_length
                newread.reference_start = read.query_alignment_start + start
                newread.reference_end = read.query_alignment_end
                newread.reference_name = chrom

                pos = alignmentstart+1
                char = 0

                newcigar = []
                clip_5 = 0
                clip_3 = 0

                changed = False
                inseq = False
                for op, length in read.cigar:
                    if op == 5:  # H
                        changed = True
                    elif op == 4:  # S
                        changed = True
                        if not inseq:
                            clip_5 = length
                        else:
                            clip_3 = length
                    else:
                        inseq = True
                        newcigar.append((op, length))

                newread.cigar = newcigar

            out.write(newread)

    close_bam(samfile)


def parse_bam(bam):
    if '.bam' in bam:
        return pysam.AlignmentFile(bam, "rb")
    elif '.sam' in bam:
        return pysam.AlignmentFile(bam, "r")


def read_head(bam):
    if '.bam' in bam:
        return pysam.AlignmentFile(bam, "rb").header
    elif '.sam' in bam:
        return pysam.AlignmentFile(bam, "r").header


def close_bam(samfile):
    samfile.close()


def check_idx(file):
    if (os.path.isfile(file+'.idx')) or (os.path.isfile(str.join('.',split('.',file)[0:-1],'fai'))) or (os.path.isfile(str.join('.',split('.',file)[0:-1],'faidx'))):
        return True
    else:
        return False


###MAIN###

if __name__ == '__main__':
    args=parseargs()
    process(args.fasta, args.bams, args.outdir, args.cluster)

#
# RemoveSoftClip.py ends here
