# Collection.py ---
#
# Filename: Collection.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Mon Jul  1 17:54:01 2019 (+0200)
#           By: Joerg Fallmann
#     Update #: 220
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
#import os, sys, inspect
# # realpath() will make your script run, even if you symlink it :)
# cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) ))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)
#
# # Use this if you want to include modules from a subfolder
# cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"lib")
# if cmd_subfolder not in sys.path:
#     sys.path.insert(0, cmd_subfolder)
#
# # Info:
# # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# # __file__ fails if the script is called in different ways on Windows.
# # __file__ fails if someone does os.chdir() before.
# # sys.argv[0] also fails, because it doesn't not always contains the path.

import glob, os, snakemake
import numpy as np
import heapq
from operator import itemgetter
from natsort import natsorted, ns
import traceback as tb
import sys
import re
import pprint
from io import StringIO
from Bio import SeqIO
import gzip
import math
import collections

# Code:All subs from here on

##############################
########Snakemake Subs########
##############################
#def sources(config):
#    try:
#        ret = list()
#        for x,y in config["SOURCE"].items():
#            if x in list(config["SAMPLES"]):
#                ret.append(str(y))
#        return ret
#
#    except Exception as err:
#        exc_type, exc_value, exc_tb = sys.exc_info()
#        tbe = tb.TracebackException(
#            exc_type, exc_value, exc_tb,
#        )
#        with open('error','a') as h:
#            print(''.join(tbe.format()), file=h)

def sources(config):
    try:
        ret = list()
        for key in config["SOURCE"] and config["SAMPLES"]:
            ret.append(str(key))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def samples(config):
    try:
        ret = list()
        for x,y in config["SOURCE"].items():
            k = find_innermost_value_from_dict(config["SAMPLES"][x])
            for l in k:
                ret.append(os.path.join(str(x),str(l)))
            print('sret ',ret)
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def sampleslong(config):
    try:
        ret = list()
        for x,y in config["SOURCE"].items():
            for s in config["SAMPLES"][x]:
                k = list_all_values_of_dict(config["SAMPLES"][x][s])
                for v in k:
                    ret.append(os.path.join(str(x),str(s),str(v)))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def samplesonly(config):        # THIS IS NOT ADVISED, SAMPLES INDEPENDENT OF SOURCE!
    try:
        ret = list()
        for x,y in config["SOURCE"].items():
            for s in config["SAMPLES"][x]:
                for n in config["SAMPLES"][x][s]:
                    ret.append(str(n))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def get_placeholder(config):
    try:
        ret = list()
        if 'PH' in (config):
            for x in config['PH']:
                ret.append(str(x))
        else:
            ret.append('_')
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def genomepath(s, config):
    try:
        s = os.path.basename(str(s))
        for k,v in config["SAMPLES"].items():
            for g,l in v.items():
                if s in l:
                    for a,b in config["SOURCE"].items():
                        for c,d in b.items():
                            if g == c:
                                for x, y in config["GENOME"].items():
                                    if str(d) == str(y):
                                        return os.path.join(str(x),str(y))
                                    elif str(d) == str(x):
                                        return os.path.join(str(x),str(y))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def genome(s, config):
    try:
        s = os.path.basename(str(s))
        for k,v in config["SAMPLES"].items():
            for g,l in v.items():
                if s in l:
                    for a,b in config["SOURCE"].items():
                        for c,d in b.items():
                            if g == c:
                                for x, y in config["GENOME"].items():
                                    if str(d) == str(y):
                                        return str(y)
                                    elif str(d) == str(x):
                                        return str(x)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def fullgenomepath(sa, config):
    try:
        ret=list()
        for s in sa:
            l = config["GENOME"][s]
            ret.append(os.path.join(str(s),str(l)))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def genomename(s, config):
    try:
        s = os.path.basename(str(s))
        for k,v in config["SAMPLES"].items():
            for g,l in v.items():
                if s in l:
                    for x, y in config["GENOME"].items():
                        if g == y:
                            return str(x)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def pathstogenomes(samples, config):
    try:
        ret = list()
        for s in samples:
            s = os.path.basename(s)
            for k,v in config["SAMPLES"].items():
                for g,l in v.items():
                    if s in l:
                        for x, y in config["GENOME"].items():
                            if g == y:
                                ret.append(os.path.join(str(x),str(y)))
        return sorted(list(set(ret)))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def mapping_params(sample, runstate, config):
    try:
        s = os.path.basename(str(sample))
        t = genome(s,config)
        mp = list()
        if runstate is None:
            runstate = str(os.path.split(runstate_from_sample(str(sample)))[-1])
        for k,v in config["MAPPING"].items():
            for g,p in v[runstate].items():
                if str(g) == str(t):
                    mp.append()
        return mp
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def count_params(sample, config):
    try:
        s = os.path.basename(str(sample))
        t = genome(s,config)
        for k,v in config["COUNT"].items():
            for g,p in v.items():
                if g == t:
                    return str(p)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def trnascan_params(s, runstate, config):
    try:
        for k,v in config["TRNASCAN"].items():
            for g,p in v[runstate].items():
                if g == s:
                    return str(p)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def index(w, t, config):
    try:
        gen = genome(str(w),config)
        return expand("{ref}/{gen}.{name}_all_withoutPseudo_cluster.idx",ref=REFERENCE,gen=gen, name=NAME[t])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def rmempty(check):
    try:
        ret = list()
        for f in check:
            if os.path.isfile(f):
                ret.append(f)
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def source_from_sample(sample):
    try:
        ret = str(os.path.join(*os.path.split(str(sample))[0:-1]))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def sample_from_path(path):
    try:
        ret = str(os.path.join(os.path.split(str(path))[-1]))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def anno_from_file(s, config):
    try:
        s = str(s)
        for k,v in config["ANNOTATION"].items():
            for g,l in v.items():
                if s == g:
                    return str(l)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def runstate_from_sample(sample):
    try:
        runstate = list()
        for s in samples:
            s = os.path.basename(s)
            for k,v in config["SAMPLES"].items():
                for g,l in list_all_values_of_dict(v):
                        for s in l:
                            for x, y in config["GENOME"].items():
                                if g == y:
                                    ret.append(os.path.join(str(x),str(y)))
        return ret
#        return sorted(list(set(ret)))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


##############################
#########Python Subs##########
##############################
def list_all_values_of_dict(dictionary):
    try:
        if isinstance(dictionary, dict):
            for values in dictionary.values():
                if isinstance(values, dict):
                    yield from list_all_values_of_dict(values)
                else:
                    yield values
        else:
            yield dictionary

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def find_all_values_on_key(key, dictionary):
    try:
        if isinstance(dictionary, dict):
            for k, v in dictionary.items():
                if k == key:
                    yield v
                if isinstance(v, dict):
                    yield from find_all_values_on_key(key, v)
        else:
            yield dictionary
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def find_innermost_value_from_dict(dictionary):
    try:
        if isinstance(dictionary, dict):
            for k, v in dictionary.items():
                if isinstance(v, dict):
                     return(find_innermost_value_from_dict(v))
                else:
                    return v
        else:
            return dictionary
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def removekey(d, key):
    try:
        r = dict(d)
        del r[key]
        return r
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def getlowest_list(a, n):
    try:
        if n > len(a) - 1:
            b = len(a) - 1
        else:
            b = n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[:n])
        else:
            return list(None for i in range(n))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def gethighest_list(a, n):
    try:
        if len(a)-n < 0:
            b = len(a)-1
        else:
            b = len(a)-n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[-n:])
        else:
            return list(None for i in range(n))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def getlowest_dict(a, n):
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nsmallest(b,a.items(), key=itemgetter(1)))
        else:
            return dict({i:None for i in range(n)})
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def gethighest_dict(a, n):
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nlargest(b,a.items(), key=itemgetter(1)))
        else:
            return dict({i:None for i in range(n)})
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def toarray(file, ulim):
    try:
        x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
        return x
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def convertcol(entry):
    try:
        if isinvalid(entry):
#       if entry is None or entry == 'NA' or entry == 'nan' or entry is np.nan:
            return np.nan
        else:
            return float(entry)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def eprint(*args, **kwargs):
    try:
        with open('error', 'a') as e:
            print(*args, file=e, **kwargs)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def printlog(msg):
    try:
        with open ('log','a') as l:
            print(str(msg), file = l )
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def parse_annotation_bed(bed, annotated=None):
    try:
        anno = defaultdict(list)
        if os.path.isfile(os.path.abspath(bed)):
            if '.gz' in bed:
                f = gzip.open(bed,'rt')
            else:
                f = open(bed,'rt')
        else:
            f = bed
        for line in f:
            entries = line.rstrip().split('\t')
            goi = entries[3]
            if annotated:
                start = int(entries[10])
                end   = int(entries[11])-1
            else:
                start = int(entries[1])
                end   = int(entries[2])-1
            anno[str(goi)].append('-'.join([str(start),str(end)]))
        return anno
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def readConstraintsFromBed(bed, linewise=None):
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            start = int(entries[1])+1
            end = entries[2]
            goi = entries[3]
            if linewise:
                cons['lw'].append('-'.join([str(start),str(end)]))
            else:
                cons[str(goi)].append('-'.join([str(start),str(end)]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def readPairedConstraintsFromBed(bed, linewise=None):
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            if entries[1] > -1 and entries[8] > -1:
                start_one = int(entries[1])+1
                end_one = entries[2]
                goi = entries[3]
                start_two = int(entries[8])+1
                end_two = entries[9]
                if linewise:
                    cons['lw'].append('-'.join([str(start_one),str(end_one),str(start_two),str(end_two)]))
                else:
                    cons[str(goi)].append('-'.join([str(start_one),str(end_one),str(start_two),str(end_two)]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def readConstraintsFromCSV(csv, linewise=None):
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = split(',',line.rstrip())
            if linewise:
                cons['def'].append('-'.join([str(entries[1]),str(entries[2])]))
            else:
                cons[entries[3]].append('-'.join([str(entries[1]),str(entries[2])]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def readConstraintsFromGeneric(generic, linewise=None):
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = re.split(r'[ ,|;"]+', line.rstrip())
            if len(entries > 2):
                if linewise:
                    cons['lw'].append('-'.join([str(entries[1]),str(entries[2])]))
                else:
                    cons[entries[0]].append('-'.join([str(entries[1]),str(entries[2])]))
            else:
                if linewise:
                    cons['lw'].append('-'.join([str(entries[1]),str(entries[2])]))
                else:
                    cons['generic'].append('-'.join(str(entries[1:2])))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def isvalid(x=None):
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return False
            else:
                return True
        else:
            return False
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def isinvalid(x=None):
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def makeoutdir(outdir):
    try:
        if not os.path.isabs(outdir):
            outdir =  os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        return outdir
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def parseseq(sequence):
    try:
        if (isinstance(sequence, StringIO)):
            seq = sequence

        elif ( isinstance(sequence, str) and sequence == 'random' ):
            rand = "\n".join(createrandseq(length, gc, number, alphabet))
            seq = StringIO(rand)
            o = gzip.open('Random.fa.gz','wb')
            o.write(bytes(rand,encoding='UTF-8'))
            o.close()

        elif (isinstance(sequence, str) and os.path.isfile(sequence)):
            if '.gz' in sequence :
                seq = gzip.open(sequence,'rt')
            else:
                seq = open(sequence,'rt')
        else:
            header = ">Seq1:default:nochrom:(.)"
            s = sequence
            seq = StringIO("{header}\n{s}".format(header=header, s=s))

        return seq
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def plot_data(fa, raw, consu, consp, const, xs, cons, saveas, outdir):
    try:
        anime = []
        #define xs for constraint line
        consl = []
        for x in const:
            consl.append(1.25)
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)

        ax2 = ax1.twiny()
    #   line, = ax.plot([], [], lw=2)
        plt.title("Blue-- = Unconstraint, Green-. = Unpaired, Red = Paired, Gray = Constraint",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
    #   plt.xticks(range(0,len(fa.seq)+1),(' '+fa.seq),size='small')
    #add lines to plot
        ax1.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax2.set_xlim(ax1.get_xlim())
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq), ha="right")
    #   ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
        ax2.set_xticks(range(1,len(fa.seq)+1))
        ax2.set_xticklabels(range(1,len(fa.seq)+1), rotation=45, ha="right")
    # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        ax2.tick_params(axis='both', which='major', labelsize=5)
        ax2.tick_params(axis='both', which='minor', labelsize=3)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('StruCons_'+goi+'_'+cons+'.'+saveas)
        plt.close()
    #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
    #   return anime
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def plot_temp(fa, raw, temp, xs, saveas, outdir):
    try:
        anime = []
        #define xs for constraint line
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)
        plt.title("Blue-- = "+temp+" degree",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
    #add lines to plot
        ax1.plot(xs, raw, 'b-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq))
    # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('TempCons_'+goi+'_'+temp+'.'+saveas)
        plt.close()
    #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
    #   return anime
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def calc_gibbs(fc):
    try:
        return fc.pf()[1]
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def get_bppm(tmp, start, end):
    try:
        bppm = []
        for item in tmp:
            for i in range(int(start),int(end)+1):
                try:
                    if item[i] > 0.0:
                        bppm.append(str.join('\t',[str(tmp.index(item)), str(i), str(item[i])]))
                except Exception as err:
                    exc_type, exc_value, exc_tb = sys.exc_info()
                    tbe = tb.TracebackException(
                        exc_type, exc_value, exc_tb,
                        )
                    with open('error','a') as h:
                        print(''.join(tbe.format()), file=h)
        return bppm
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def calc_bpp(bppm):
    bpp = 0.0
    try:
        for entry in bppm:
            base, mate, prob = map(float,entry.split('\t'))
            bpp += prob
        return bpp

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def calc_nrg(bpp):
    try:
        #set kT for nrg2prob and vice versa calcs
        kT = 0.61632077549999997

        nrg = 0.0;
        if bpp > 0.0:
            nrg = -1 * kT * math.log(bpp)
        return nrg
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def print_region_up(data, seqlength=None, region=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    try:
        if data:
            ups=''
            x = int(region)
            for i in range(int(seqlength)):
                if isinvalid(data[i][x]):
                    data[i][x] = np.nan
                else:
                    data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+str(data[i][x])+"\n"
            return ups
        else:
            eprint('No up data to print')
            return ups

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def print_up(data=None, seqlength=None, region=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    try:
        if data:
            ups=''
            for i in range(int(seqlength)):
                for x in range(1,region+1):
                    if isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
            return ups
        else:
            eprint('No up data to print')
            return ups
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


def up_to_array(data=None, region=None, seqlength=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data[165553:165588])
    try:
        if data:
            entries=[]
            if not seqlength:
                seqlength = len(data)
            if not region:
                region = slice(1,len(data[0]))
            for i in range(seqlength):
                entries.append([])
                for e in range(len(data[i])):
                    if isinvalid(data[i][e]):
                        data[i][e] = np.nan
                    else:
                        data[i][e] = round(data[i][e],8)
                entries[i].extend(data[i][region])
            return np.array(entries)
        else:
            eprint('No up data to print')
            return np.array()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def npprint(a, o=None):#, format_string ='{0:.2f}'):
    try:
        out = ''
        it = np.nditer(a, flags=['f_index'])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index+1,it[0])+"\n"
            it.iternext()
        if o:
            o.write(bytes(out,encoding='UTF-8'))
        else:
            print(out)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def printdiff(a, o=None):
    try:
        np.savetxt(o, a, delimiter='\t')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def read_precalc_plfold(data, name, seq):
    try:
        for i in range(len(seq)):
            data.append([])
            data[i] = []
        with gzip.open(name,'rt') as o:
            for line in o:
                cells = line.rstrip().split('\t')
                data[int(cells[0])-1].append([])
                data[int(cells[0])-1][0] = None
                for a in range(1,len(cells)):
                    data[int(cells[0])-1].append([])
                    data[int(cells[0])-1][a] = float(cells[a])
        return data
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)
    return 1

def pl_to_array(name, ulim):
    try:
        printlog(name)
        return np.array(np.loadtxt(name, usecols=ulim, unpack=True, delimiter='\t'))
#        data = []
#        with gzip.open(name,'rt') as o:
#            for line in o:
#                cells = line.rstrip().split('\t')
#                data.append(float(cells[ulim]))
#        return data
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def idfromfa(id):
    goi, chrom, strand = [None, None, None]
    try:
        goi, chrom = id.split(':')[::2]
        strand = str(id.split(':')[3].split('(')[1][0])
    except:
        eprint('Fasta header is not in expected format, you will loose information on strand and chromosome')
        goi = id
        chrom, strand = ['na','na']

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        sys.exit('Could not assign any value from fasta header, please check your fasta files')

def cluster2trna(seqs):
    try:
        translater = collections.OrderedDict()
        translater['cluster'] = collections.OrderedDict()
        translater['tRNA'] = collections.OrderedDict()

        for fa in SeqIO.parse(seqs,'fasta'):
            head = str(fa.id).upper()           # cluster1:chr19.tRNA5-LysCTT(+) cluster2:NC_007091.3.tRNA25-ArgTCT
            cluster, info = head.split(':')
            chrom, trna = (info.split('.')[0], info.split('.')[-1].split('(')[0])
            strand = 'u'
            if '(+)' in info or '(-)' in info:
                strand = re.sub('[()]', '_', info.split('.')[-1]).split('_')[1]

            if cluster in translater['cluster']:
                translater['cluster'][cluster].append(trna)
            else:
                translater['cluster'][cluster] = list()
                translater['cluster'][cluster].append(trna)
            if chrom in translater['tRNA']:
                if strand in translater['tRNA'][chrom]:
                    translater['tRNA'][chrom][strand].append(trna)
                else:
                    translater['tRNA'][chrom][strand] = list()
                    translater['tRNA'][chrom][strand].append(trna)
            else:
                translater['tRNA'][chrom] = collections.OrderedDict()
                translater['tRNA'][chrom][strand] = list()
                translater['tRNA'][chrom][strand].append(trna)

        return translater

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)

def check_ref(reference):
    try:
        if os.path.exists(os.path.abspath(reference)):
            return reference
        elif os.path.exists(os.path.abspath(reference+'.gz')):
            return reference+'.gz'

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        with open('error','a') as h:
            print(''.join(tbe.format()), file=h)


#
# Collection.py ends here
