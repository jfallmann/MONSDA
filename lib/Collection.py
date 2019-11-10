# Collection.py ---
#
# Filename: Collection.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Sun Nov 10 12:02:48 2019 (-0300)
#           By: Joerg Fallmann
#     Update #: 785
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
import inspect
import subprocess
import collections
from collections import defaultdict, OrderedDict
import six

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

log = setup_logger(name='Collection', log_file='LOGS/Snakemake_Collection.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')

#Class
class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))

# Code:All subs from here on
##############################
########Snakemake Subs########
##############################
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
        log.error(''.join(tbe.format()))

def samples(config):
    try:
        ret = list()
        for x,y in config["SOURCE"].items():
            k = find_innermost_value_from_dict(config["SAMPLES"][x])
            for l in k:
                ret.append(os.path.join(str(x),str(l)))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def sampleslong(config):
    try:
        ret = list()
        for x,y in config["SOURCE"].items():
            for s in config["SAMPLES"][x]:
                k = list_all_values_of_dict(config["SAMPLES"][x][s])
                for v in k:
                    if isinstance(v, list):
                        for z in v:
                            ret.append(os.path.join(str(x),str(s),str(z)))
                    else:
                        ret.append(os.path.join(str(x),str(s),str(v)))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))


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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))


def genomepath(s, config):
    try:
        sa = os.path.basename(str(s))
        cond= s.split(os.sep)[-2]
        sk = find_key_for_value(sa, config["SAMPLES"])
        for skey in sk:
            klist = value_extract(skey, config["SOURCE"])
            for k in klist:
                for x, y in config["GENOME"].items():
                    if str(k) == str(y) or str(k) == str(x):
                        return os.path.join(str(x),str(y))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def genome(s, config):
    try:
        logid='genome: '
        sa = os.path.basename(str(s))
        sp = source_from_sample(str(s)).split(os.sep)[0]
        cond= s.split(os.sep)[-2]
        sk = find_key_for_value(sa, config["SAMPLES"])
        for skey in sk:
            klist = value_extract(skey, config["SOURCE"])
            for k in klist:
                if str(k) == sp:
                    log.debug(logid+'k is sp')
                    for x, y in config["GENOME"].items():
                        log.debug(logid+str([k, x, y]))
                        if str(k) == str(x):
                            return str(y)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def namefromfile(s, config):
    try:
        sa = os.path.basename(str(s))
        cond= s.split(os.sep)[-2]
        sk = find_key_for_value(sa, config["SAMPLES"])
        for skey in sk:
            klist = value_extract(skey, config["NAME"])
            for k in klist:
                    if str(skey) == str(cond):
                        return str(k)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def create_subworkflow(config, subwork, conditions):
    try:
        logid = 'create_subworkflow: '
        toollist = list()
        configs = list()
        for condition in conditions:
            try:
                env = str(subDict(config[subwork],condition)['ENV'])
            except:
                log.warning('Key ENV not found for '+subwork+' this can be intentional')
                env = ''
            try:
                exe = str(subDict(config[subwork],condition)['BIN'])
            except:
                log.warning('Key BIN not found for '+subwork+' this can be intentional')
                exe = ''
            src, treat, setup = condition
            log.debug(logid+str([env,exe,src,treat,setup]))
            tempconf = NestedDefaultDict()
            try:
                for key in ['REFERENCE', 'BINS','MAXTHREADS']:
                    tempconf[key] = config[key]
            except KeyError:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(
                    exc_type, exc_value, exc_tb,
                )
                log.error(''.join(tbe.format()))
            try:
                for key in ['GENOME', 'NAME']:
                    tempconf[key][src] = config[key][src]
                for key in ['SOURCE', 'SAMPLES', 'SEQUENCING', subwork]:
                    tempconf[key][src][treat][setup] = config[key][src][treat][setup]
            except KeyError:
                exc_type, exc_value, exc_tb = sys.exc_info()
                tbe = tb.TracebackException(
                    exc_type, exc_value, exc_tb,
                )
                log.error(''.join(tbe.format()))
            tempconf[subwork+'ENV'] = env
            tempconf[subwork+'BIN'] = exe
            toollist.append([env,exe])
            configs.append(tempconf)
        return toollist, configs

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def namefrompath(p, config):
    try:
        p = os.path.dirname(p).split(os.sep)
        klist = getFromDict(config["NAME"],p)
        for k in klist:
            return str(k)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

#def tool_params(sample, runstate, config, subconf):
#    try:
#        logid='tool_params: '
#        t = genome(sample,config)
#        mp = list()
#        if runstate is None:
#            runstate = runstate_from_sample([sample], config)
#        x = source_from_sample(sample).split(os.sep)
#        for k in getFromDict(config[subconf],x):
#            log.debug(logid+str(k))
#            y = find_key_for_value(k,config[subconf])
#            for r in runstate:
#                if r in sample.split(os.sep) and r in [z for z in y]:
#                    mp.extend(k)
#        log.debug(logid+str(mp))
#        return mp
#    except Exception as err:
#        exc_type, exc_value, exc_tb = sys.exc_info()
#        tbe = tb.TracebackException(
#            exc_type, exc_value, exc_tb,
#        )
#        log.error(''.join(tbe.format()))

#def index_params(indexpath, config, subconf):
#    try:
#        s = indexpath.split(os.sep)
#        mp = list()
#        for k in getFromDict(config[subconf],s):
#            mp.extend(k)
#        return mp
#    except Exception as err:
#        exc_type, exc_value, exc_tb = sys.exc_info()
#        tbe = tb.TracebackException(
#            exc_type, exc_value, exc_tb,
#        )
#        log.error(''.join(tbe.format()))

def tool_params(sample, runstate, config, subconf):
    try:
        logid='tool_params: '
        log.debug(logid+'Samples: '+str(sample))
        t = genome(sample,config)
        mp = OrderedDict()
        if runstate is None:
            runstate = runstate_from_sample([sample], config)[0]
        x = source_from_sample(sample).split(os.sep)
        log.debug(logid+str([sample,runstate,subconf,t,x]))
        if runstate not in x:
            x.append(runstate)
        mp = subDict(config[subconf],x)
        log.debug(logid+'DONE: '+str(mp))
        return mp
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def env_bin_from_config(samples, config, subconf):
    try:
        s = samples[0].split(os.sep)[:-1]
        mb,me = [None,None]
        for k in getFromDict(config[subconf],s):
            mb, me = k['BIN'], k['ENV']
        return mb,me
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def env_bin_from_config2(samples, config, subconf):
    try:
        logid='env_bin_from_config2'
        for s in samples:
            log.debug(logid+': '+s)
            log.debug(str(config[subconf]))
            for k in getFromDict(config[subconf],conditiononly(s,config)):
                mb = k['BIN']
                me = k['ENV']
        log.debug([str(mb),str(me)])
        return mb, me
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def index(w, t, config):
    try:
        gen = genome(str(w),config)
        return expand("{ref}/{gen}.{name}_all_withoutPseudo_cluster.idx",ref=REFERENCE,gen=gen, name=NAME[t])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def source_from_sample(sample):
    try:
        ret = str(os.path.join(*os.path.split(str(sample))[0:-1]))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def sample_from_path(path):
    try:
        ret = str(os.path.join(os.path.split(str(path))[-1]))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def anno_from_file(sample, config, step):
    try:
        logid = 'anno_from_file: '
        p = os.path.dirname(genomepath(sample, config))
        s = source_from_sample(sample)
        ret = os.path.join(config["REFERENCE"],p,subDict(config["ANNOTATE"],s)[step])
        log.debug(logid+str(ret))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def anno_from_source(source, config, step):
    try:
        logid = 'anno_from_source: '
        s = source.split(os.sep)[0:-1]
        p = s[0]
        samp = source.split(os.sep)[-1]
        log.debug(logid+str(s))
        runstate = runstate_from_sample([samp], config)[0]
        lst = list()
        lst.extend(s)
        lst.append(runstate)
        log.debug(logid+str(lst))
        ret = os.path.join(config["REFERENCE"],p,subDict(config["ANNOTATE"],lst)[step])
        log.debug(logid+str(ret))
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def runstate_from_sample(sample,config):
    try:
        logid = 'runstate_from_sample: '
        ret = list()
        for s in sample:
            s = os.path.basename(s)
            for k,v in config["SAMPLES"].items():
                for f in find_key_for_value(s,v):
                    log.debug(logid+f)
                    ret.append(f)
        return ret
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def samplecond(sample,config):
    try:
        logid = 'samplecond: '
        ret = list()
        paired = False
        for s in sample:
            check = os.path.dirname(s).split(os.sep)
            log.debug(logid+str(check))
            for r in runstate_from_sample([s],config):
                tmplist = check
                tmplist.append(r)
                if getFromDict(config['SEQUENCING'],tmplist) is 'paired':
                    paired = True
                if paired:
                    s=re.sub(r'_[r|R|\A\Z][1|2]','',s)
                ret.append(os.path.join("{p}".format(p=os.path.dirname(s)),"{c}".format(c=r),os.path.basename(s)))
        log.debug(logid+str([sample,ret]))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def conditiononly(sample,config):
    try:
        logid = 'conditiononly: '
        ret = list()
        paired = False
        check = os.path.dirname(sample).split(os.sep)
        log.debug(logid+str(check))
        for r in runstate_from_sample([sample],config):
            ret.extend(check)
            ret.append(r)
        log.debug(logid+str([sample,ret]))
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def checkpaired(sample,config):
    try:
        logid = 'checkpaired: '
        ret = list()
        paired = False
        for s in sample:
            check = os.path.dirname(s).split(os.sep)
            log.debug(logid+str(check))
            for r in runstate_from_sample([s],config):
                tmplist = check
                tmplist.append(r)
                if 'paired' in getFromDict(config['SEQUENCING'],tmplist):
                    paired = True
        log.debug(logid+str(paired))
        return paired

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def checkclip(sample,config):
    try:
        logid = 'checkclip: '
        ret = list()
        clip = ''
        for s in sample:
            check = os.path.dirname(s).split(os.sep)
            log.debug(logid+str(check))
            r = runstate_from_sample([s],config)
            tmplist = check
            tmplist.extend(r)
            log.debug(logid+str(tmplist))
            if 'CLIP' in subDict(config['PEAKS'],tmplist):
                clip = subDict(config['PEAKS'],tmplist)['CLIP']
            else:
                log.debug(logid+'Key CLIP not found in config')
        log.debug(logid+str(clip))
        return str(clip)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def aggregate_input(wildcards):
    return expand("post/{sample}/{i}.txt",
           sample=wildcards.sample,
           i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)

##############################
#########Python Subs##########
##############################
def dict_inst(d):
    try:
        loginfo='dict_inst: '
        if isinstance(d,dict) or isinstance(d,OrderedDict) or isinstance(d,defaultdict):
            return True
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def getFromDict(dataDict, mapList):
    logid = 'getFromDict: '
    log.debug(logid+str(mapList))
    ret=list()
    for k in mapList:
        dataDict = dataDict[k]
    ret.append(dataDict)
    log.debug(logid+str(ret))
    return ret

def subDict(dataDict, mapList):
    ret=dict()
    for k in mapList:
        dataDict = dataDict[k]
    return dataDict

def nested_set(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value

def merge_dicts(d,u):
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


def list_all_keys_of_dict(dictionary):
    try:
        if dict_inst(dictionary):
            for key in dictionary.keys():
                if dict_inst(key):
                    yield from list_all_keys_of_dict(key)
                else:
                    yield key
        else:
            yield dictionary

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def list_all_values_of_dict(dictionary):
    try:
        if dict_inst(dictionary):
            for values in dictionary.values():
                if dict_inst(values):
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
        log.error(''.join(tbe.format()))

def find_all_values_on_key(key, dictionary):
    try:
        if dict_inst(dictionary):
            for k, v in dictionary.items():
                if dict_inst(v):
                    yield from find_all_values_on_key(key, v)
                elif k == key:
                    yield v

        else:
            return dictionary
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def find_key_for_value(val, dictionary):
    try:
        logid='find_key_for_value: '
        log.debug(logid+str(val))
        if dict_inst(dictionary):
            for k, v in dictionary.items():
                if dict_inst(v):
                    log.debug(logid+'item'+str(v))
                    yield from find_key_for_value(val, v)
                elif v == val or val in v:
                    yield k
        else:
            return dictionary
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def value_extract(key, var):
    try:
        logid='value_extract: '
        log.debug(logid+str(var))
        if hasattr(var,'items'):
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
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

def find_innermost_value_from_dict(dictionary):
    try:
        if dict_inst(dictionary):
            for k, v in dictionary.items():
                if dict_inst(v):
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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def toarray(file, ulim):
    try:
        x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
        return x
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))


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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def calc_gibbs(fc):
    try:
        return fc.pf()[1]
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
                    log.error(''.join(tbe.format()))
        return bppm
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))


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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))


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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def printdiff(a, o=None):
    try:
        np.savetxt(o, a, delimiter='\t')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))
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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

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
        log.error(''.join(tbe.format()))

def runjob(jobtorun):
    return subprocess.run(jobtorun, shell=True, universal_newlines=True, capture_output=True)  # python >= 3.7

#
# Collection.py ends here
