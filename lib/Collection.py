# Collection.py ---
#
# Filename: Collection.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Fri Apr 17 20:53:57 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 1783
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
# # sys.argv[0] also fails, because it doesn't not always contain the path.

import glob, os, snakemake
import numpy as np
import heapq
import itertools
from itertools import combinations
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
import shlex
import collections
from collections import defaultdict, OrderedDict
import six
import logging

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename)
    if any(x in scriptname for x in ['Snakemake','Configurator']):
        log = logging.getLogger(scriptname)
    else:
        log = logging.getLogger('SubSnake')
        lvl = log.level
        if (log.hasHandlers()):
            log.handlers.clear()
        handler = logging.FileHandler('LOGS/RunSnakemake.log',mode='a')
        handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s',datefmt='%m-%d %H:%M'))
        log.setLevel('DEBUG')
        log.addHandler(handler)  # streamlog

    log.debug('FILENAME: '+str(scriptname))

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()))

#Class
class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))

    def merge(self, *args):
        self = merge_dicts(self,*args)

# Code:All subs from here on
##############################
########Snakemake Subs########
##############################
def check_run(func):
    def func_wrapper(*args, **kwargs):
        logid=scriptname+'.Collection_func_wrapper: '
        try:
            return func(*args, **kwargs)

        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
            )
            log.error(logid+''.join(tbe.format()))
    return func_wrapper

@check_run
def sources(config):
    logid = scriptname+'.Collection_sources: '
    ret = list()
    search =  [x[0] for x in keysets_from_dict(config["SOURCE"]) if x[0] != 'last']
    if len(getFromDict(config['SAMPLES'],search)) > 0:
        ret.extend(search)
    log.debug(logid+str(ret))
    return ret

@check_run
def get_samples(config):
    logid = scriptname+'.Collection_get_samples: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid+'SAMPLES_LONG: '+str(SAMPLES))
    check = [os.path.join('FASTQ',str(x)+'*.fastq.gz') for x in SAMPLES]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        paired = checkpaired([SAMPLES[i]],config)
        log.debug(logid+'PAIRED: '+str(paired))
        log.debug(logid+'SEARCHING: '+s)
        f = glob.glob(s)
        log.debug(logid+'SAMPLECHECK: '+str(f))
        if f:
            f = list(set([str.join(os.sep,s.split(os.sep)[1:]) for s in f]))
            if paired == 'paired':
                RETSAMPLES.extend(list(set([os.path.join(os.path.dirname(s),re.sub(r'_r1|_R1|_r2|_R2|.fastq.gz','',os.path.basename(s))) for s in f])))
                log.debug(logid+'PAIREDSAMPLES: '+str(f))
            else:
                RETSAMPLES.extend([x.replace('.fastq.gz','') for x in f])
    log.debug(logid+'SAMPLETEST: '+str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid+'No samples found, please check config file')
        sys.exit()

    log.info(logid+'SAMPLES: '+str(RETSAMPLES))
    return RETSAMPLES

@check_run
def download_samples(config):
    logid = scriptname+'.Collection_download_samples: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid+'DOWNLOAD_SAMPLES_LONG: '+str(SAMPLES))
    return SAMPLES

@check_run
def get_conditions(samples, config):
    logid = scriptname+'.Collection_conditions: '
    ret = list()
    for k in keysets_from_dict(config['SOURCE']):
        ret.append(k)
    log.debug(logid+str(ret))
    return ret

@check_run
def get_samples_from_dir(id, condition, setting, config):
    logid = scriptname+'.Collection_get_samples_from_dir: '
    pat = os.path.abspath(os.path.join('FASTQ',id, condition, '*.fastq.gz'))
    log.info(logid+str(pat))
    ret = natsorted(glob.glob(pat), key=lambda y: y.lower())
    log.debug(logid+str(ret))
    if len(ret) > 0:
        seqtype = getFromDict(config, ['SEQUENCING', id, condition, setting])
        for x in seqtype:
            if 'unpaired' not in x:
                ret = list(set([re.sub(r'_r1|_R1|_r2|_R2|.fastq.gz','',os.path.basename(s)) for s in ret]))
                renamelist = [re.sub(r'_r\d', lambda pat: pat.group(1).upper(), s) for s in ret]
                for i in range(len(renamelist)):
                    if renamelist[i] != ret[i]:
                        os.rename(ret[i],renamelist[i])
            else:
                ret = list(set([re.sub(r'.fastq.gz','',os.path.basename(s)) for s in ret]))
        return list(set(ret))
    else:
        return list()

@check_run
def sampleslong(config):
    logid = scriptname+'.Collection_sampleslong: '
    ret = list()
    tosearch = list()
    for k in keysets_from_dict(config['SAMPLES']):
        tosearch.append(k)
    log.debug(logid+'keys: '+str(tosearch))
    for search in tosearch:
        for x in list(set(getFromDict(config['SAMPLES'],search)[0])):
            ret.append(os.path.join(str.join(os.sep,search[:-1]),x))
    ret= list(set(ret))
    log.debug(logid+str(ret))
    return ret

@check_run
def samplesonly(config):        # THIS IS NOT ADVISED, SAMPLES INDEPENDENT OF SOURCE!
    ret = list()
    for x,y in config["SOURCE"].items():
        for s in config["SAMPLES"][x]:
            for n in config["SAMPLES"][x][s]:
                ret.append(str(n))
    return ret

@check_run
def get_placeholder(config):
    ret = list()
    if 'PH' in (config):
        for x in config['PH']:
            ret.append(str(x))
    else:
        ret.append('_')
    return ret

@check_run
def genomepath(s, config):
    logid=scriptname+'.Collection_genomepath: '
    sa = os.path.basename(str(s))
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config["SAMPLES"])
    log.debug(logid+'GENOMEPATH: '+str([sa,cond,sk]))
    for skey in sk:
        klist = value_extract(skey, config["SOURCE"])
        for k in klist:
            for x, y in config["GENOME"].items():
                log.debug(logid+'GENOMEPATH: '+str([x,y,k]))
                if str(k) == str(y) or str(k) == str(x):
                    return os.path.join(str(x),str(y))

@check_run
def genome(s, config):
    logid=scriptname+'.Collection_genome: '
    sa = os.path.basename(str(s))
    sp = source_from_sample(str(s),config)
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config['SAMPLES'])
    for skey in sk:
        klist = value_extract(skey, config['SOURCE'])
        for k in klist:
            for x, y in config['GENOME'].items():
                log.debug(logid+str([k, x, y]))
                if str(k) == str(x):
                    return str(y)

@check_run
def fullgenomepath(sa, config):
    ret=list()
    for s in sa:
        l = config["GENOME"][s]
        ret.append(os.path.join(str(s),str(l)))
    return ret

@check_run
def genomename(s, config):
    s = os.path.basename(str(s))
    for k,v in config["SAMPLES"].items():
        for g,l in v.items():
            if s in l:
                for x, y in config["GENOME"].items():
                    if g == y:
                        return str(x)

@check_run
def namefromfile(s, config):
    if 'NAME' not in config:
        return ''
    else:
        sa = os.path.basename(str(s))
        cond= s.split(os.sep)[-2]
        sk = find_key_for_value(sa, config["SAMPLES"])
        for skey in sk:
            klist = value_extract(skey, config["NAME"])
            for k in klist:
                if str(skey) == str(cond):
                    return str(k)

@check_run
def create_subworkflow(config, subwork, conditions, stage=''):
    logid = scriptname+'.Collection_create_subworkflow: '
    log.debug(logid+str([config, subwork, conditions, stage]))
    toollist = list()
    configs = list()
    for condition in conditions:
        try:
            env = str(subDict(config[subwork],condition)[stage+'ENV'])
        except:
            if subwork not in ['DE', 'DEU', 'DAS']:
                log.warning('Key ENV not found for '+subwork+' this can be intentional')
            env = ''
        try:
            exe = str(subDict(config[subwork],condition)[stage+'BIN'])
        except:
            if subwork not in ['DE', 'DEU', 'DAS', 'RAW']:
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
            matchinggenome=config['SOURCE'][src][treat][setup]
            tempconf['GENOME'][matchinggenome] = config['GENOME'][matchinggenome]
            for key in ['NAME', 'SOURCE', 'SAMPLES', 'SEQUENCING', subwork]:
                tempconf[key][src][treat][setup] = config[key][src][treat][setup]
            if any([subwork == x for x in ['DE','DEU','DAS','COUNTING']]):
                if subwork == 'COUNTING':
                    tempconf['COUNTING']['FEATURES'] = config['COUNTING']['FEATURES']
                if subwork == 'DAS':
                    tempconf['MAPPING'][src][treat][setup] = config['MAPPING'][src][treat][setup]
                if 'COMPARABLE' in config[subwork]:
                    tempconf[subwork]['COMPARABLE'] = config[subwork]['COMPARABLE']
                if 'TOOLS' in config[subwork]:
                    tempconf[subwork]['TOOLS'] = config[subwork]['TOOLS']

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

    log.debug(logid+str([toollist,configs]))

    return toollist, configs

@check_run
def namefrompath(p, config):
    p = os.path.dirname(p).split(os.sep)
    klist = getFromDict(config["NAME"],p) if 'NAME' in config else list('')
    for k in klist:
        return str(k)

@check_run
def pathstogenomes(samples, config):
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

@check_run
def tool_params(sample, runstate, config, subconf):
    logid=scriptname+'.Collection_tool_params: '
    log.debug(logid+'Samples: '+str(sample))
    t = genome(sample,config)
    mp = OrderedDict()
    x = sample.split(os.sep)[:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    log.debug(logid+str([sample,runstate,subconf,t,x]))
    mp = subDict(config[subconf],x)
    log.debug(logid+'DONE: '+str(mp))
    return mp

@check_run
def get_reps(samples,config,analysis):
    logid=scriptname+'.Collection_get_reps: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(sample.split(os.sep)[2:-1]))
        partconf = subDict(config[analysis],sample.split(os.sep)[2:-1])
        log.debug(logid+'CONF: '+str(partconf))
        ret['reps'].append(sample)
        wcfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        idx = partconf['REPLICATES'].index(wcfile)
        ret['pairs'].append(checkpaired_rep([str.join(os.sep,sample.split(os.sep)[2:])],config))
        ret['conds'].append(partconf['GROUPS'][idx])
        if 'TYPES' in partconf and len(partconf['TYPES']) >= idx:
                ret['types'].append(partconf['TYPES'][idx])
        else:
            ret['types'].append('std')

    rets = '-r '+str.join(',',ret['reps'])
    rets += ' -c '+str.join(',',ret['conds'])
    rets += ' -t '+str.join(',',ret['types'])
    rets += ' --paired '+str.join(',',ret['pairs']) if 'pairs' in ret else ''

    log.debug(logid+'RETURN: '+str(rets))
    return rets

@check_run
def get_diego_samples(samples,config,analysis):
    logid=scriptname+'.Collection_get_diego_samples: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(sample.split(os.sep)[2:-1]))
        partconf = subDict(config[analysis],sample.split(os.sep)[2:-1])
        log.debug(logid+'CONF: '+str(partconf))
        wcfile = str.join('-',sample.split(os.sep)[-3:]).replace('_mapped_sorted_unique.counts','')
        ret[wcfile].append(sample)

    log.debug(logid+'RETURN: '+str(ret))

    slist = ''
    for key,val in ret.items():
        slist +=  str(key)+'\t'
        slist += '|'.join(val)
        slist += os.linesep

        log.debug(logid+'RETURN: '+str(slist))
    return slist

@check_run
def get_diego_groups(samples,config,analysis):
    logid=scriptname+'.Collection_get_diego_groups: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(sample.split(os.sep)[2:-1]))
        partconf = subDict(config[analysis],sample.split(os.sep)[2:-1])
        log.debug(logid+'CONF: '+str(partconf))
        wcfile = str.join('-',sample.split(os.sep)[-3:]).replace('_mapped_sorted_unique.counts','')
        checkfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        #wcfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        idx = partconf['REPLICATES'].index(checkfile)
        cond = partconf['GROUPS'][idx]
        ret[cond].append(wcfile)

    slist = ''
    for key,val in ret.items():
        slist +=  str(key)+'\t'
        slist += '|'.join(val)
        slist += os.linesep
    log.debug(logid+'RETURN: '+str(slist))
    return slist

@check_run
def env_bin_from_config(samples, config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config: '
    s = samples[0].split(os.sep)[:-1]
    mb,me = [None,None]
    for k in getFromDict(config[subconf],s):
        mb, me = k['BIN'], k['ENV']
    return mb,me

@check_run
def env_bin_from_config2(samples, config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config2: '
    for s in samples:
        log.debug(logid+'S: '+str(s))
        log.debug(logid+'C: '+str(conditiononly(s,config)))
        check = conditiononly(s,config)

        for k in getFromDict(config[subconf],check):
            if 'BIN' in k:
                mb = k['BIN']
            else:
                mb = ''
            if 'ENV' in k:
                me = k['ENV']
            else:
                me = ''
        log.debug(logid+str([str(mb),str(me)]))
    return mb, me

@check_run
def env_bin_from_config3(config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config3: '
    envkey = subconf+"ENV"
    binkey = subconf+"BIN"
    me = config[envkey]
    mb = config[binkey]
    log.debug(logid+str([str(mb),str(me)]))
    return mb, me

@check_run
def rmempty(check):
    ret = list()
    for f in check:
        if os.path.isfile(f):
            ret.append(f)
    return ret

@check_run
def source_from_sample(sample, config):
    logid=scriptname+'.Collection_source_from_sample: '
    s = os.path.dirname(str(sample))
    cond= s.split(os.sep)
    log.debug(logid+str([s,cond]))
    ret = getFromDict(config["SOURCE"],cond)[0]
    return ret

@check_run
def sample_from_path(path):
    ret = str(os.path.join(os.path.split(str(path))[-1]))
    return ret

@check_run
def anno_from_file(sample, config, step):
    logid = scriptname+'.Collection_anno_from_file: '
    p = os.path.dirname(genomepath(sample, config))
    s = source_from_sample(sample,config)
    ret = os.path.join(config["REFERENCE"],p,subDict(config["ANNOTATE"],s)[step])
    log.debug(logid+str(ret))
    return ret

@check_run
def anno_from_source(source, config, step):
    logid = scriptname+'.Collection_anno_from_source: '
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

@check_run
def runstate_from_sample(sample,config):
    logid = scriptname+'.Collection_runstate_from_sample: '
    ret = list()
    for s in sample:
        n = s.split(os.sep)[-1]
        s = os.path.dirname(s)
        if len(s.split(os.sep)) > 2:
            s = str.join(os.sep,s.split(os.sep)[-3:])
        log.debug(logid+'SAMPLE: '+s)
        c = getFromDict(config["SAMPLES"],s.split(os.sep))[0]
        if dict_inst(c):
            for k,v in c.items():
                log.debug(logid+'k,v: '+str([str(k),str(v)]))
                if n in v:
                    if k not in ret:
                        ret.append(k)
        else:
            if n in c:
                k = s.split(os.sep)[-1]
                ret.append(k)
    log.debug(logid+'RETURN: '+str(ret))
    return ret

@check_run
def samplecond(sample,config):
    logid = scriptname+'.Collection_samplecond: '
    ret = list()
    for s in sample:
        s = s.replace('.fastq.gz','')
        log.debug(logid+'SAMPLE: '+str(s))
        check = os.path.dirname(s).split(os.sep)
        log.debug(logid+'CHECK: '+str(check))
        for r in runstate_from_sample([s],config):
            tmplist = check
            tmplist.append(r)
            log.debug(logid+'TMPLIST: '+str(tmplist))
            if getFromDict(config['SEQUENCING'],tmplist) is 'paired':
                s=re.sub(r'_[r|R|\A\Z][1|2]','',s)
            ret.append(os.path.join("{p}".format(p=os.path.dirname(s)),"{c}".format(c=r),os.path.basename(s)))
    log.debug(logid+'RETURN: '+str(ret))
    return ret

@check_run
def conditiononly(sample,config):
    logid = scriptname+'.Collection_conditiononly: '
    ret = list()
    paired = False
    check = os.path.dirname(sample).split(os.sep)
    log.debug(logid+str(check))
    for r in runstate_from_sample([sample],config):
        log.debug(logid+str(r))
        if len(ret) < 3:
            ret.extend(check)
            if r not in ret:
                if len(ret) < 3:
                    ret.append(r)
    log.debug(logid+str(ret))
    return ret

@check_run
def checkpaired(sample,config):
    logid = scriptname+'.Collection_checkpaired: '
    ret = list()
    paired = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(s))
        check = os.path.dirname(s).split(os.sep)
        tmplist = check
        p = getFromDict(config['SEQUENCING'],tmplist)[0]
        log.debug(logid+'P: '+str(p))
        for r in runstate_from_sample([s],config):
            if r in p:
                tmplist.append(r)
                paired = getFromDict(config['SEQUENCING'],tmplist)[0].split(',')[0]
                tmplist = tmplist[:2]
    log.debug(logid+'PAIRED: '+str(paired))
    return paired

@check_run
def checkpaired_rep(sample,config):
    logid = scriptname+'.Collection_checkpaired_rep: '
    log.debug(logid+'SAMPLE: '+str(sample))
    ret = list()
    for s in sample:
        check = os.path.dirname(s).split(os.sep)
        tmplist = check
        log.debug(logid+'S: '+str(tmplist))
        p = getFromDict(config['SEQUENCING'],tmplist)[0]
        ret.append(p)
    log.debug(logid+'PAIRED: '+str(ret))
    return str.join(',',ret)

@check_run
def checkstranded(sample,config):
    logid = scriptname+'.Collection_checkstranded: '
    ret = list()
    stranded = ''
    for s in sample:
        check = os.path.dirname(s).split(os.sep)
        tmplist = check
        p = getFromDict(config['SEQUENCING'],tmplist)[0]
        log.debug(logid+'P: '+str(p))
        for r in runstate_from_sample([s],config):
            if r in p:
                tmplist.append(r)
                stranded = getFromDict(config['SEQUENCING'],tmplist)[0].split(',')[1] if len(getFromDict(config['SEQUENCING'],tmplist)[0].split(',')) > 1 else ''
                tmplist = tmplist[:2]
    log.debug(logid+'STRANDEDNESS: '+str(stranded))
    return stranded

@check_run
def post_checkpaired(sample,config):
    logid = scriptname+'.Collection_checkpaired: '
    ret = list()
    paired = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(sample))
        check = os.path.dirname(s).split(os.sep)
        tmplist = check
        log.debug(logid+'TMP: '+str(tmplist))
        p = getFromDict(config['SEQUENCING'],tmplist)[0]
        log.debug(logid+'P: '+str(p))
        #if not dict_inst(p):
        #paired = p[0] if 'paired' in p or 'unpaired' in p or 'singlecell' in p else ''
        log.debug(logid+'P: '+str(p))
        for r in runstate_from_sample([s],config):
            log.debug(logid+'R: '+str(r))
            if r in p:
                tmplist.append(r)
                paired = getFromDict(config['SEQUENCING'],tmplist)[0].split(',')[0]
                tmplist = tmplist[:2]
    log.debug(logid+'PAIRED: '+str(paired))
    return paired

@check_run
def checkclip(sample,config):
    logid = scriptname+'.Collection_checkclip: '
    ret = list()
    clip = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(s))
        check = os.path.dirname(s).split(os.sep)
        r = runstate_from_sample([s],config)
        log.debug(logid+'RUNSTATE: '+str(r))
        tmplist = check
        if r not in tmplist:
            tmplist.extend(r)
        log.debug(logid+str(tmplist))
        log.debug(logid+'TMP: '+str(tmplist))
        check = getFromDict(config['PEAKS'],tmplist)[0]
        log.debug(logid+'CHECK: '+str(check))
        if 'CLIP' in check:
            clip = check['CLIP']
        else:
            log.debug(logid+'Key CLIP not found in config')
    log.debug(logid+'CLIP is: '+str(clip))
    return str(clip)

@check_run
def check_tool_params(sample, runstate, config, subconf, idx):
    try:
        par = tool_params(sample, runstate ,config, subconf)['OPTIONS'][idx]
        if par is not '':
            return par
        elif subconf == 'MAPPING':
            return 'std'
        else:
            return ''
    except:
        if subconf == 'MAPPING':
            return 'std'
        else:
            return ''

@check_run
def comparable_as_string(config, subwork):
    logid=scriptname+'.comparable_as_string: '
    check = config[subwork].get('COMPARABLE')
    if check:
        log.info(logid+'determine comparables in '+subwork)
        complist  = []
        compdict=config[subwork]['COMPARABLE']
        for key in compdict:
            for value in compdict[key]:
                complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr
    else:
        log.info(logid+'no comparables found in '+subwork+'. Compare All vs. All.')
        groups_by_condition = list(yield_from_dict("GROUPS",config))
        flattened = set(val for sublist in groups_by_condition for val in sublist)
        combined=list(combinations(flattened,2))
        complist=[]
        for key, value in combined:
            complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr

@check_run
def comparable_as_string2(config, subwork):
    logid=scriptname+'.comparable_as_string2: '
    check = config[subwork].get('COMPARABLE')
    if check:
        log.info(logid+'determine comparables in '+subwork)
        complist  = []
        compdict=config[subwork]['COMPARABLE']
        for contrast in compdict:
            As = ""
            Bs = ""
            for condition in compdict[contrast][0]:
                As = (As + "+" + condition).strip("+")
            for condition in compdict[contrast][1]:
                Bs = (Bs + "+" + condition).strip("+")
            complist.append(f"{contrast}:{As}-vs-{Bs}")
        compstr = ','.join(complist)
        return compstr
    else:
        log.info(logid+'no comparables found in '+subwork+'. Compare All vs. All.')
        groups_by_condition = list(yield_from_dict("GROUPS",config))
        flattened = set(val for sublist in groups_by_condition for val in sublist)
        combined=list(combinations(flattened,2))
        complist=[]
        for key, value in combined:
            complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr

##############################
#########Python Subs##########
##############################
@check_run
def dict_inst(d):
    logid=scriptname+'.Collection_dict_inst: '
    if isinstance(d,dict) or isinstance(d,OrderedDict) or isinstance(d,defaultdict) or isinstance(d,NestedDefaultDict):
        return True

@check_run
def getFromDict(dataDict, mapList):
    logid = scriptname+'.Collection_getFromDict: '
    log.debug(logid+'MAPLIST: '+str(mapList))
    ret = dataDict
    for k in mapList:
        if k in dataDict:
            log.debug(logid+'k: '+str(k))
            dataDict = dataDict[k]
        else:
            return list([])
    log.debug(logid+'MIDRET: '+str(ret))
    if ret != dataDict:
        log.debug(logid+'RET: '+str(dataDict))
        return list([dataDict])
    else:
        log.debug(logid+'RET: '+str(list([])))
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
    logid = scriptname+'.Collection_subDict: '
    log.debug(logid+str(mapList))
    ret = dataDict
    for k in mapList:
        log.debug(logid+'k: '+str(k))
        if k in ret:
            ret = ret[k]
    return ret

@check_run
def nested_set(dic, keys, value):
    for key in keys[:-1]:
        dic = dic.setdefault(key, {})
    dic[keys[-1]] = value

@check_run
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

@check_run
def keysets_from_dict(dictionary,original=None):  # Only works for equal depth keysets, needs tweaking for other use cases
    logid = scriptname+'.Collection_keysets_from_dict: '

    keylist = list()
    if dict_inst(dictionary):
        for k,v in keys_from_dict(dictionary).items():
            keylist.append(v)
        log.debug(logid+'kl:'+str(keylist))
        combis = list(itertools.product(*keylist))
        log.debug(logid+'cs:'+str(combis))
        ret = list()
        for combi in combis:
            log.debug(logid+'combi: '+str(combi))
            if len(getFromDict(dictionary,combi)) >= 1:
                log.debug(logid+'found: '+str(combi))
                ret.append(combi)
            else:
                continue
        return ret
    else:
        return keylist

@check_run
def keys_from_dict(dictionary,first=True,lvl=0,save=None):
    logid = scriptname+'.Collection_keys_from_dict: '

    if first:
        first = False
        end = depth(dictionary)
        save = defaultdict(list)
        log.debug(logid+'END:'+str(end))

    if dict_inst(dictionary):
        log.debug(logid+'dictDEPTH: '+str(depth(dictionary)))
        log.debug(logid+'dictLEVEL: '+str(lvl))
        for k,v in dictionary.items():
            save[lvl].append(k)
            log.debug(logid+str(save))
            if dict_inst(v):
                save = keys_from_dict(v,first,lvl+1,save)
            else:
                continue
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
    logid = scriptname+'.Collection_list_all_keys_of_dict: '
    for key, value in dictionary.items():
        if type(value) is dict:
            yield key
            yield from recursive_items(value)
        else:
            yield key

@check_run
def list_all_values_of_dict(dictionary):
    if dict_inst(dictionary):
        for key, value in dictionary.items():
            if dict_inst(value):
                #yield (key, value)
                yield from list_all_values_of_dict(value)
            else:
                yield (key, value)
                #yield ('last','value')
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
    logid=scriptname+'.Collection_find_key_for_value: '
    log.debug(logid+'VAL: '+str(val)+' Dict: '+str(dictionary))
    if dict_inst(dictionary):
        for k, v in dictionary.items():
            log.debug(logid+'DICT: '+str(k)+' ; '+str(v))
            if dict_inst(v):
                log.debug(logid+'item'+str(v))
                yield from find_key_for_value(val, v)
            elif v == val or val in v:
                yield k
    else:
        return dictionary

@check_run
def value_extract(key, var):
    logid=scriptname+'.Collection_value_extract: '
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

@check_run
def find_innermost_value_from_dict(dictionary):
    logid=scriptname+'.Collection_find_innermost_value_from_dict: '
    if dict_inst(dictionary):
        for k, v in dictionary.items():
            if dict_inst(v):
                 return(find_innermost_value_from_dict(v))
            else:
                return v
    else:
        return dictionary
    return ret

@check_run
def removekey(d, key):
    logid=scriptname+'.Collection_removekey: '
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
    if len(a)-n < 0:
        b = len(a)-1
    else:
        b = len(a)-n
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
        return dict(heapq.nsmallest(b,a.items(), key=itemgetter(1)))
    else:
        return dict({i:None for i in range(n)})

@check_run
def gethighest_dict(a, n):
    if n > len(a):
        b = len(a)
    else:
        b = n
    if len(a) > 0:
        return dict(heapq.nlargest(b,a.items(), key=itemgetter(1)))
    else:
        return dict({i:None for i in range(n)})

@check_run
def toarray(file, ulim):
    x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
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
        if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
            return False
        else:
            return True
    else:
        return False

@check_run
def isinvalid(x=None):

    if x:
        if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
            return True
        else:
            return False
    else:
        return True

@check_run
def makeoutdir(outdir):

    if not os.path.isabs(outdir):
        outdir =  os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

@check_run
def parseseq(sequence):
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

@check_run
def npprint(a, o=None):#, format_string ='{0:.2f}'):
    out = ''
    it = np.nditer(a, flags=['f_index'])
    while not it.finished:
        out += "%d\t%0.7f" % (it.index+1,it[0])+"\n"
        it.iternext()
    if o:
        o.write(bytes(out,encoding='UTF-8'))
    else:
        print(out)

@check_run
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

@check_run
def cluster2trna(seqs):
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

@check_run
def check_ref(reference):
    if os.path.exists(os.path.abspath(reference)):
        return reference
    elif os.path.exists(os.path.abspath(reference+'.gz')):
        return reference+'.gz'

def multi_replace(repl, text):
    print('MULTI: '+str(repl)+str(text))
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, repl.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)

#
# Collection.py ends here
