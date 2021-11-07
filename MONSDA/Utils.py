# Utils.py ---
#
# Filename: Utils.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Feb  4 18:01:07 2021 (+0100)
#           By: Joerg Fallmann
#     Update #: 2888
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
# import os, sys, inspect
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
# # sys.argv[0] also fails, because it doesn't not always contain the path.

import os
import sys
import re
import glob
import shutil
import json
import numpy as np
import heapq
import itertools
from operator import itemgetter
from natsort import natsorted
import traceback as tb
from io import StringIO
from Bio import SeqIO
import gzip
import inspect
import subprocess
import collections
from collections import defaultdict, OrderedDict
import six
import logging
import hashlib
from snakemake import load_configfile
import functools
import datetime

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


#  NestedDefaultDict
class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))

    def merge(self, *args):
        self = merge_dicts(self, *args)


# Wrapper
def check_run(func):
    @functools.wraps(func)
    def func_wrapper(*args, **kwargs):
        logid = scriptname + ".Collection_func_wrapper: "
        try:
            return func(*args, **kwargs)

        except Exception:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type,
                exc_value,
                exc_tb,
            )
            log.error(logid + "".join(tbe.format()))

    return func_wrapper


@check_run
def rmempty(check):
    ret = list()
    for f in check:
        if os.path.isfile(f):
            ret.append(f)
    return ret


@check_run
def ns_check_version(v, r):
    logid = scriptname + ".MONSDA_check_version: "

    if parse_version(v) < parse_version(check):
        log.debug(logid + check)
        return True
    else:
        return shutil.which("MONSDA")


##############################
#########Python Subs##########
##############################
@check_run
def dict_inst(d):
    logid = scriptname + ".Collection_dict_inst: "
    if (
        isinstance(d, dict)
        or isinstance(d, OrderedDict)
        or isinstance(d, defaultdict)
        or isinstance(d, NestedDefaultDict)
    ):
        return True


@check_run
def getFromDict(dataDict, mapList):
    logid = scriptname + ".Collection_getFromDict: "
    log.debug(logid + "MAPLIST: " + str(mapList) + "\tDict: " + str(dataDict))
    ret = dataDict
    for k in mapList:
        log.debug(logid + "k: " + str(k))
        if dataDict.get(k):
            dataDict = dataDict[k]
            log.debug(logid + "subdict: " + str(dataDict))
        else:
            return list([])
    if ret != dataDict:
        log.debug(logid + "RET: " + str(dataDict))
        return list([dataDict])
    else:
        log.debug(logid + "RET: " + str(list([])))
        return list([])


@check_run
def yield_from_dict(key, dictionary):
    for k, v in dictionary.items():
        if k == key:
            yield v
        elif dict_inst(v):
            for result in yield_from_dict(key, v):
                yield result
        elif isinstance(v, list):
            for d in v:
                if dict_inst(d):
                    for result in yield_from_dict(key, d):
                        yield result


@check_run
def subDict(dataDict, mapList):
    logid = scriptname + ".Collection_subDict: "
    log.debug(logid + str(mapList))
    ret = dataDict
    for k in mapList:
        log.debug(logid + "k: " + str(k))
        if k in ret:
            ret = ret[k]
            if not dict_inst(ret):
                ret = {k: ret}
        else:
            log.debug(logid + "No k in dict")
            return dict()
    return ret


@check_run
def subSetDict(dataDict, mapList):
    logid = scriptname + ".Collection_subSetDict: "
    log.debug(logid + str(mapList))
    parse = subDict(dataDict, mapList)
    ret = {}
    nested_set(ret, mapList, parse)
    log.debug(logid + str(ret))
    return ret


@check_run
def nested_set(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value


@check_run
def merge_dicts(d, u):
    # https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    # python 3.8+ compatibility
    try:
        collectionsAbc = collections.abc
    except:
        collectionsAbc = collections

    for k, v in six.iteritems(u):
        dv = d.get(k, {})
        if not isinstance(dv, collectionsAbc.Mapping):
            d[k] = v
        elif isinstance(v, collectionsAbc.Mapping):
            d[k] = merge_dicts(dv, v)
        else:
            d[k] = v
    return d


@check_run
def keysets_from_dict(
    dictionary, search=None, original=None
):  # Only works for equal depth keysets, needs tweaking for other use cases
    logid = scriptname + ".Collection_keysets_from_dict: "

    keylist = list()
    if dict_inst(dictionary):
        for k, v in keys_from_dict(dictionary, search).items():
            keylist.append(v)
        log.debug(logid + "kl:" + str(keylist))
        combis = list()
        for i in range(1, len(keylist) + 1):
            subkeylist = keylist[0:i]
            combis.extend(list(itertools.product(*subkeylist)))
        log.debug(logid + "cs:" + str(combis))
        ret = list()
        for combi in combis:
            # if len(getFromDict(dictionary, combi)) >= 1:
            check = subDict(dictionary, combi)
            log.debug(logid + "checking: " + str(check))
            if isvalid(check):
                if (isinstance(check, dict) and check.get("SAMPLES")) or isinstance(
                    check, str
                ):
                    log.debug(logid + "found: " + str(combi))
                    ret.append(combi)
        return ret
    else:
        return keylist


@check_run
def keys_from_dict(dictionary, search=None, save=None, first=True, lvl=0):
    logid = scriptname + ".Collection_keys_from_dict: "

    if first:
        first = False
        end = depth(dictionary)
        save = defaultdict(list)
        log.debug(logid + "TOTALDEPTH:" + str(end))

    if dict_inst(dictionary):
        log.debug(logid + "dictDEPTH: " + str(depth(dictionary)))
        log.debug(logid + "dictLEVEL: " + str(lvl))
        for k, v in dictionary.items():
            if not search or (search and k != search):
                save[lvl].append(k)
                log.debug(logid + "TMPSAVE: " + str(save))
                if dict_inst(v):
                    save = keys_from_dict(v, search, save, first, lvl + 1)
                else:
                    continue
            else:
                log.debug(logid + "Found search: " + str(save))
                return save
        return save
    else:
        return save


@check_run
def depth(d):
    if dict_inst(d):
        return 1 + (max(map(depth, d.values())) if d else 0)
    return 0


@check_run
def list_all_keys_of_dict(dictionary):
    logid = scriptname + ".Collection_list_all_keys_of_dict: "
    for key, value in dictionary.items():
        if type(value) is dict:
            yield key
            yield from list_all_keys_of_dict(value)
        else:
            yield key


@check_run
def list_all_values_of_dict(dictionary):
    if dict_inst(dictionary):
        for key, value in dictionary.items():
            if dict_inst(value):
                # yield (key, value)
                yield from list_all_values_of_dict(value)
            else:
                yield (key, value)
                # yield ('last','value')
    else:
        yield dictionary


@check_run
def find_all_values_on_key(key, dictionary):
    if dict_inst(dictionary):
        for k, v in dictionary.items():
            if dict_inst(v):
                yield from find_all_values_on_key(key, v)
            elif k == key:
                yield v

    else:
        return dictionary


@check_run
def find_key_for_value(val, dictionary):
    logid = scriptname + ".Collection_find_key_for_value: "
    log.debug(logid + "VAL: " + str(val) + " Dict: " + str(dictionary))
    if dict_inst(dictionary):
        for k, v in dictionary.items():
            log.debug(logid + "DICT: " + str(k) + " ; " + str(v))
            if dict_inst(v):
                log.debug(logid + "item" + str(v))
                yield from find_key_for_value(val, v)
            elif v == val or val in v:
                yield k
    else:
        return dictionary


@check_run
def value_extract(key, var):
    logid = scriptname + ".Collection_value_extract: "
    log.debug(logid + str(var))
    if hasattr(var, "items"):
        for k, v in var.items():
            if k == key:
                yield v
            if dict_inst(v):
                for result in value_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in value_extract(key, d):
                        yield result


@check_run
def find_innermost_value_from_dict(dictionary):
    logid = scriptname + ".Collection_find_innermost_value_from_dict: "
    if dict_inst(dictionary):
        for k, v in dictionary.items():
            if dict_inst(v):
                return find_innermost_value_from_dict(v)
            else:
                return v
    else:
        return dictionary


@check_run
def removekey(d, key):
    logid = scriptname + ".Collection_removekey: "
    r = dict(d)
    del r[key]
    return r


@check_run
def getlowest_list(a, n):
    if n > len(a) - 1:
        b = len(a) - 1
    else:
        b = n
    if len(a) > 0 and n > 0:
        return list(np.partition(a, b)[:n])
    else:
        return list(None for i in range(n))


@check_run
def gethighest_list(a, n):
    if len(a) - n < 0:
        b = len(a) - 1
    else:
        b = len(a) - n
    if len(a) > 0 and n > 0:
        return list(np.partition(a, b)[-n:])
    else:
        return list(None for i in range(n))


@check_run
def getlowest_dict(a, n):
    if n > len(a):
        b = len(a)
    else:
        b = n
    if len(a) > 0:
        return dict(heapq.nsmallest(b, a.items(), key=itemgetter(1)))
    else:
        return dict({i: None for i in range(n)})


@check_run
def gethighest_dict(a, n):
    if n > len(a):
        b = len(a)
    else:
        b = n
    if len(a) > 0:
        return dict(heapq.nlargest(b, a.items(), key=itemgetter(1)))
    else:
        return dict({i: None for i in range(n)})


@check_run
def toarray(file, ulim):
    x = np.loadtxt(
        str(file),
        usecols=(ulim),
        delimiter="\t",
        unpack=True,
        converters={ulim: lambda s: convertcol(s.decode("utf-8"))},
    )
    return x


@check_run
def convertcol(entry):
    if isinvalid(entry):
        return np.nan
    else:
        return float(entry)


@check_run
def isvalid(x=None):
    if x:
        if x in ("None", "nan", "none", "NA", "NAN") or x is None or x is np.nan:
            return False
        else:
            return True
    else:
        return False


@check_run
def isinvalid(x=None):

    if x:
        if x in ("None", "nan", "none", "NA", "NAN") or x is None or x is np.nan:
            return True
        else:
            return False
    else:
        return True


@check_run
def makeoutdir(outdir):

    if not os.path.isabs(outdir):
        outdir = os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir


@check_run
def parseseq(sequence):
    if isinstance(sequence, StringIO):
        seq = sequence

    elif isinstance(sequence, str) and os.path.isfile(sequence):
        if ".gz" in sequence:
            seq = gzip.open(sequence, "rt")
        else:
            seq = open(sequence, "rt")
    else:
        header = ">Seq1:default:nochrom:(.)"
        s = sequence
        seq = StringIO("{header}\n{s}".format(header=header, s=s))

    return seq


@check_run
def npprint(a, o=None):  # , format_string ='{0:.2f}'):
    out = ""
    it = np.nditer(a, flags=["f_index"])
    while not it.finished:
        out += "%d\t%0.7f" % (it.index + 1, it[0]) + "\n"
        it.iternext()
    if o:
        o.write(bytes(out, encoding="UTF-8"))
    else:
        print(out)


@check_run
def idfromfa(id):
    goi, chrom, strand = [None, None, None]
    try:
        goi, chrom = id.split(":")[::2]
        strand = str(id.split(":")[3].split("(")[1][0])
    except:
        print(
            "Fasta header is not in expected format, you will loose information on strand and chromosome"
        )
        goi = id
        chrom, strand = ["na", "na"]

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        sys.exit(
            "Could not assign any value from fasta header, please check your fasta files"
        )


@check_run
def cluster2trna(seqs):
    translater = collections.OrderedDict()
    translater["cluster"] = collections.OrderedDict()
    translater["tRNA"] = collections.OrderedDict()

    for fa in SeqIO.parse(seqs, "fasta"):
        head = str(
            fa.id
        ).upper()  # cluster1:chr19.tRNA5-LysCTT(+) cluster2:NC_007091.3.tRNA25-ArgTCT
        cluster, info = head.split(":")
        chrom, trna = (info.split(".")[0], info.split(".")[-1].split("(")[0])
        strand = "u"
        if "(+)" in info or "(-)" in info:
            strand = re.sub("[()]", "_", info.split(".")[-1]).split("_")[1]

        if cluster in translater["cluster"]:
            translater["cluster"][cluster].append(trna)
        else:
            translater["cluster"][cluster] = list()
            translater["cluster"][cluster].append(trna)
        if chrom in translater["tRNA"]:
            if strand in translater["tRNA"][chrom]:
                translater["tRNA"][chrom][strand].append(trna)
            else:
                translater["tRNA"][chrom][strand] = list()
                translater["tRNA"][chrom][strand].append(trna)
        else:
            translater["tRNA"][chrom] = collections.OrderedDict()
            translater["tRNA"][chrom][strand] = list()
            translater["tRNA"][chrom][strand].append(trna)

    return translater


@check_run
def check_ref(reference):
    if os.path.exists(os.path.abspath(reference)):
        return reference
    elif os.path.exists(os.path.abspath(reference + ".gz")):
        return reference + ".gz"


@check_run
def multi_replace(repl, text):
    print("MULTI: " + str(repl) + str(text))
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, repl.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start() : mo.end()]], text)


@check_run
def makelogdir(logdir):
    if not os.path.isabs(logdir):
        logdir = os.path.abspath(logdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    return logdir


@check_run
def get_dict_hash(d):
    logid = scriptname + ".get_dict_hash: "
    log.debug(logid + "INPUT DICT: " + str(d))
    ret = str(hashlib.sha256(bytes(str(sorted(d.items())), "utf-8")).hexdigest())
    log.debug(logid + "HASH: " + ret)
    return ret


@check_run
def add_to_innermost_key_by_list(addto, toadd, keylist):
    logid = scriptname + ".add_to_innermost_key_by_list: "
    log.debug(logid + str(addto) + ", " + str(toadd))

    tconf = {}
    for i in range(
        len(keylist)
    ):  # need to add options as last element to dict of unknown depth
        tconf[keylist[i]] = {}
        tconf = tconf[keylist[i]]
        if i == len(keylist) - 1:
            tconf = toadd

    addto.update(tconf)
    return addto


#
# Utils.py ends here
