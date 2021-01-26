# Collection.py ---
#
# Filename: Collection.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Tue Sep 18 15:39:06 2018 (+0200)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Jan 26 14:58:10 2021 (+0100)
#           By: Joerg Fallmann
#     Update #: 2773
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
# TODO for arbitrary depth:
# runstate_from_sample, sampleslong, samplecond,
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

import glob
import os
import snakemake
import json
import numpy as np
import heapq
import itertools
from itertools import combinations, product
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
import hashlib
from snakemake import load_configfile
import functools
import json
import textwrap
# from string import join, split


try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace('.py','')
    if any(x in scriptname for x in ['RunSnakemake','Configurator']):
        log = logging.getLogger(scriptname)
    else:
        log = logging.getLogger('snakemake')

    lvl = log.level if log.level else 'INFO'
    for handler in log.handlers[:]:
        handler.close()
        log.removeHandler(handler)

    handler = logging.FileHandler('LOGS/RunSnakemake.log', mode='a')
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    log.setLevel(lvl)

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    print(''.join(tbe.format()), file=sys.stderr)


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
    @functools.wraps(func)
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
def get_samples(config):
    logid = scriptname+'.Collection_get_samples: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid+'SAMPLES_LONG: '+str(SAMPLES))
    check = [os.path.join('FASTQ', str(x).replace('.fastq.gz', '')+'*.fastq.gz') for x in SAMPLES]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        log.debug(logid+'SEARCHING: '+s)
        paired = checkpaired([SAMPLES[i]], config)
        log.debug(logid+'PAIRED: '+str(paired))
        f = glob.glob(s)
        log.debug(logid+'SAMPLECHECK: '+str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            if 'paired' not in paired:
                RETSAMPLES.extend(list(set([os.path.join(os.path.dirname(s), re.sub(r'_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz', '', os.path.basename(s))) for s in f])))
                log.debug(logid+'PAIREDSAMPLES: '+str(f))
            else:
                RETSAMPLES.extend([x.replace('.fastq.gz', '') for x in f])
        else:
            log.debug(logid+'SAMPLECHECK: '+str(f))
    log.debug(logid+'SAMPLETEST: '+str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid+'No samples found, please check config file')
        sys.exit()

    log.debug(logid+'SAMPLES: '+str(RETSAMPLES))
    return RETSAMPLES


@check_run
def get_samples_postprocess(config, subwork):
    logid = scriptname+'.Collection_get_samples_postprocess: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config) if len(getFromDict(config[subwork], samplecond(x, config))) > 0 ] #only those samples where postprocessing steps are defined for in config
    log.debug(logid+'SAMPLES_LONG: '+str(SAMPLES))
    check = [os.path.join('FASTQ', str(x).replace('.fastq.gz','')+'*.fastq.gz') for x in SAMPLES]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        paired = checkpaired([SAMPLES[i]], config)
        log.debug(logid+'PAIRED: '+str(paired))
        log.debug(logid+'SEARCHING: '+s)
        f = glob.glob(s)
        log.debug(logid+'SAMPLECHECK: '+str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            if 'paired' in paired:
                RETSAMPLES.extend(list(set([os.path.join(os.path.dirname(s), re.sub(r'_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz','', os.path.basename(s))) for s in f])))
                log.debug(logid+'PAIREDSAMPLES: '+str(f))
            else:
                RETSAMPLES.extend([x.replace('.fastq.gz', '') for x in f])
    log.debug(logid+'SAMPLETEST: '+str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid+'No samples found for '+str(subwork)+', please check config file')
        sys.exit()

    log.debug(logid+'SAMPLES: '+str(RETSAMPLES))
    return RETSAMPLES

@check_run
def download_samples(config):
    logid = scriptname+'.Collection_download_samples: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid+'DOWNLOAD_SAMPLES_LONG: '+str(SAMPLES))
    return SAMPLES

@check_run
def basecall_samples(config):
    logid = scriptname+'.Collection_basecall_samples: '
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    log.debug(logid+'SAMPLES_LONG: '+str(SAMPLES))
    check = [os.path.join('RAW', str(x).replace('.fast5','')+'*.fast5') for x in SAMPLES]
    RETSAMPLES = list()
    for i in range(len(check)):
        s = check[i]
        log.debug(logid+'SEARCHING: '+s)
        f = glob.glob(s)
        log.debug(logid+'SAMPLECHECK: '+str(f))
        if f:
            f = list(set([str.join(os.sep, s.split(os.sep)[1:]) for s in f]))
            RETSAMPLES.extend([x.replace('.fast5','') for x in f])
    log.debug(logid+'SAMPLETEST: '+str(RETSAMPLES))
    if len(RETSAMPLES) < 1:
        log.error(logid+'No samples found, please check config file')
        sys.exit()

    log.debug(logid+'SAMPLES: '+str(RETSAMPLES))
    return RETSAMPLES


@check_run
def get_conditions(samples, config):
    logid = scriptname+'.Collection_conditions: '
    ret = list()
    for k in keysets_from_dict(config['SETTINGS'], 'SAMPLES'):
        ret.append(k)
    log.debug(logid+str(ret))
    return list(set(ret))


@check_run
def get_samples_from_dir(search, config):  # CHECK
    logid = scriptname+'.Collection_get_samples_from_dir: '
    samples = [x.replace(' ', '') for x in list(set(getFromDict(config['SETTINGS'], search)[0]['SAMPLES']))]
    for x in range(len(search), len(search)-2, -1):  # For arbitrary depth of ics we append subdirectories until samples are found, maximum of one setting additional to file path is allowed
        pat = os.sep.join(['FASTQ', os.sep.join(search[0:x]), '*.fastq.gz'])
        log.debug(logid+'REGEX: '+str(pat)+'\t'+'SAMPLES: '+str(samples))
        check = natsorted(glob.glob(pat), key=lambda y: y.lower())
        log.debug(logid+'check: '+str(check))
        if len(check) > 0:
            ret = list()
            paired = checkpaired([os.sep.join([os.sep.join(search), check[0].split(os.sep)[-1]])], config)
            if 'paired' in paired:
                log.debug(logid+'SEARCHING: '+str([os.sep.join([os.sep.join(os.path.dirname(s).split(os.sep)[1:]), re.sub(r'_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz', '', os.path.basename(s))]) for s in check]))
                ret.extend(list(set([os.sep.join([os.sep.join(os.path.dirname(s).split(os.sep)[1:]), re.sub(r'_r1.fastq.gz|_R1.fastq.gz|_r2.fastq.gz|_R2.fastq.gz|.fastq.gz', '', os.path.basename(s))]) for s in check if any(c in s for c in samples)])))
                log.debug(logid+'FOUND: '+str(ret))
                renamelist = [re.sub(r'_r\d', lambda pat: pat.group(1).upper(), s) for s in ret]
                for i in range(len(renamelist)):
                    if renamelist[i] != ret[i]:
                        log.warning('SAMPLE NAMES CONTAIN LOWER CASE r1/r2 INSTEAD OF R1/R2 FOR PAIRED END SEQUENCING, THEY WILL BE RENAMED')
                        os.rename(ret[i], renamelist[i])
            else:
                log.debug(logid+'SEARCHING: '+str([os.sep.join(s.split(os.sep)[1:]).replace('.fastq.gz', '') for s in check]))
                ret.extend([os.sep.join(s.split(os.sep)[1:]).replace('.fastq.gz', '') for s in check if any(c in s for c in samples)])

            log.debug(logid+'RETURN: '+str(ret))
            return list(set(ret))
    log.warning(logid+'NO SAMPLES FOUND')
    return list()


@check_run
def sampleslong(config):
    logid = scriptname+'.Collection_sampleslong: '
    tosearch = list()
    for k in keysets_from_dict(config['SETTINGS'], 'SAMPLES'):  # CHECK
        tosearch.append(list(k))
    log.debug(logid+'keys: '+str(tosearch))
    for search in tosearch:
        ret = get_samples_from_dir(search, config)
    log.debug(logid+str(ret))
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
def get_cutoff_as_string(config, subwork):
    logid=scriptname+'.get_cutoff: '
    cutoffs = [f"{k}:{v}" for (k, v) in config[subwork]['CUTOFF'][0].items()]
    ret = '-'.join(cutoffs)
    log.info(logid+'CUTOFFS: '+str(ret))
    return ret


@check_run
def get_summary_dirs(config):
    logid=scriptname+'.get_summary_dirs: '
    ret=list()
    for work, tools in config['WORKFLOWS'].items():
        for tool in tools:
            ret.append(f"{work}/{tool.upper()}")
    log.debug(logid+str(ret))
    return ret


@check_run
def get_summary_files(config):
    logid=scriptname+'.get_summary_files: '
    ret=list()
    for work, tools in config['WORKFLOWS'].items():
        for tool in tools:
            log.info(logid+'make summary of '+str(work)+' - '+str(tool))
            for f in glob.glob(f"{work}/{tool.upper()}/Sig*"):
            # for f in glob.glob(f"{work}/{tool.upper()}/*"):
                ret.append(f)
    log.debug(logid+str(ret))
    return ret


@check_run
def create_subworkflow(config, subwork, conditions, stage='', combination=None):
    logid = scriptname+'.Collection_create_subworkflow: '
    log.debug(logid+str([config, subwork, conditions, stage]))
    toollist = list()
    configs = list()
    for condition in conditions:
        try:
            env = str(subDict(config[subwork], condition)[stage+'ENV'])
        except:
            if 'TOOLS' not in config[subwork]:
                log.error('No tool environment found for '+subwork+'! Either key ENV or TOOLS must be set for '+str(condition)+'!')
            env = ''
        try:
            exe = str(subDict(config[subwork], condition)[stage+'BIN'])
        except:
            if 'TOOLS' not in config[subwork]:
                log.error('No tool binary found for '+subwork+'! Either key BIN or TOOLS must be set for '+str(condition)+'!')
            exe = ''

        tempconf = NestedDefaultDict()
        tempconf['SAMPLES'] = subDict(config['SETTINGS'], condition)['SAMPLES']
        if env != '' and exe != '':
            toollist.append([env, exe])
            tempconf[subwork+'ENV'] = env
            tempconf[subwork+'BIN'] = exe

        try:
            for key in ['BINS', 'MAXTHREADS']:
                tempconf[key] = config[key]
        except KeyError:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
            )
            log.error(''.join(tbe.format()))
        try:
            for key in ['SETTINGS', subwork]:
                if len(getFromDict(config[subwork], condition)) <1:
                    if any([subwork == x for x in ['QC', 'DEDUP', 'TRIMMING', 'MAPPING']]):
                        log.error(logid+'Keys '+str(condition)+' not defined for '+str(key))
                    else:
                        log.warning(logid+'Keys '+str(condition)+' not defined for '+str(key)+', will be removed from SAMPLES for this analysis')
                        toollist.append([None, None])
                        configs.append(None)
                else:
                    tempconf[key] = subSetDict(config[key], condition)
                    if key == 'SETTINGS' and config.get('DEDUP') and 'DEDUP' in config['WORKFLOWS']:
                        tempconf['RUNDEDUP'] = 'enabled'

            if 'TOOLS' in config[subwork] and env == '' and exe == '':  # env and exe overrule TOOLS
                tempconf[subwork]['TOOLS'] = config[subwork]['TOOLS']
                for k, v in config[subwork]['TOOLS'].items():
                    toollist.append([k, v])

            if any([subwork == x for x in ['PEAKS', 'DE', 'DEU', 'DAS', 'DTU', 'COUNTING']]):
                if config[subwork].get('CUTOFF'):
                    tempconf[subwork]['CUTOFF'] = config[subwork]['CUTOFF']  #else '.05'
                if subwork == 'COUNTING':
                    tempconf['COUNTING']['FEATURES'] = config['COUNTING']['FEATURES']
                if subwork == 'DAS':
                    tempconf['MAPPING'] = subsetDict(config['MAPPING'], condition)
                if 'COMPARABLE' in config[subwork]:
                    tempconf[subwork]['COMPARABLE'] = config[subwork]['COMPARABLE']

        except KeyError:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
            )
            log.error(''.join(tbe.format()))

        configs.append(tempconf)

    log.debug(logid+str([toollist, configs]))

    return toollist, configs

@check_run
def make_pre(subwork, config, samples, condition, subdir, loglevel, state='', subname=None):
    logid = scriptname+'.Collection_make_pre: '
    log.debug(logid+'WORK: '+str(subwork))

    subjobs = list()
    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
    if listoftools is None:
        log.warning(logid+'No entry fits condition '+str(condition)+' for processing step '+str(subwork))
        return None

    for i in range(0, len(listoftools)):
        toolenv, toolbin = map(str, listoftools[i])
        subconf = NestedDefaultDict()
        sconf = listofconfigs[0]
        sconf[subwork+'ENV'] = toolenv
        sconf[subwork+'BIN'] = toolbin
        subconf.update(sconf)
        subname = toolenv+'.smk'
        log.debug(logid+str(subwork)+': '+str([toolenv, subname, condition, subconf]))

        smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'header.smk'))
        with open(smkf,'r') as smk:
            for line in smk.readlines():
                line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                line = re.sub(condapath, 'conda:  "../', line)
                subjobs.append(line)
            subjobs.append('\n\n')

        if subwork == 'QC':
            subname = toolenv+'_raw.smk'

        # Add variable for combination string
        subjobs.append('\ncombo = \''+str(subname.replace('.smk', ''))+os.sep+'\'\n\nwildcard_constraints:\n\tcombo = combo,\n \trawfile = \'|\'.join(list(SAMPLES)),\n\tfile = \'|\'.join(list(samplecond(SAMPLES, config))),\n\tread = "R1|R2"\n')
        subjobs.append('\n\n')
        # Add rulethemall based on chosen workflows
        subjobs.append(''.join(rulethemall(subwork, config, loglevel, condapath, logfix, subname.replace('.smk', ''))))        # RuleThemAll for snakemake depending on chosen workflows

        smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', subname))
        with open(smkf, 'r') as smk:
            for line in smk.readlines():
                line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                line = re.sub(condapath, 'conda:  "../', line)
                subjobs.append(line)
            subjobs.append('\n\n')

        smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'footer.smk'))
        with open(smkf, 'r') as smk:
            for line in smk.readlines():
                line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                line = re.sub(condapath, 'conda:  "../', line)
                subjobs.append(line)
            subjobs.append('\n\n')

        smko = os.path.abspath(os.path.join(subdir, '_'.join(['_'.join(condition), state+subwork, toolenv, 'subsnake.smk'])))

        if os.path.exists(smko):
            os.rename(smko, smko+'.bak')
        with open(smko, 'a') as smkout:
            smkout.write(''.join(subjobs))
            smkout.write('\n\n')

        confo = os.path.abspath(os.path.join(subdir, '_'.join(['_'.join(condition), state+subwork, toolenv, 'subconfig.json'])))
        if os.path.exists(confo):
            os.rename(confo, confo+'.bak')
        with open(confo, 'a') as confout:
            json.dump(subconf, confout)

        jobs.append([smko, confo])

    log.debug(logid+'JOBS: '+str(jobs))
    return jobs

@check_run
def get_combo(wfs, config, conditions):
    logid=scriptname+'.Collection_get_combo: '
    log.debug(logid+str(wfs)+str(conditions))
    combos = NestedDefaultDict()

    if wfs is None:
        return None

    for condition in conditions:
        ret = list()
        for subwork in wfs:
            listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
            if listoftools is None:
                log.warning(logid+'No entry fits condition '+str(condition)+' for processing step '+str(subwork))
                return None

            tools = list()
            for k, v in [toolenvs for toolenvs in listoftools]:
                tools.append({subwork:k})
                log.debug(logid+'Preparing '+str(subwork)+' for condition '+str(condition)+' with Tool: '+str(tools))
            ret.append(tools)

        combos[condition] = itertools.product(*ret)

    #no debug log here or iterator will be consumed
    return combos

@check_run
def get_combo_name(combinations):
    logid=scriptname+'.Collection_get_combo_name: '
    combname = NestedDefaultDict()

    for condition in combinations:
        combname[condition]['envs'] = list()
        combname[condition]['works'] = list()
        combos = combinations[condition]
        for combi in combos:
            envs = list()
            works = list()
            for step in combi:
                for work, env in step.items():
                    envs.append(env)
                    works.append(work)
            combname[condition]['envs'].append(str.join('-', envs))
            combname[condition]['works'].append(str.join('-', works))

    log.debug(logid+'ComboName: '+str(combname))
    return combname

@check_run
def make_sub(subworkflows, config, samples, conditions, subdir, loglevel, subname=None, combinations=None):
    logid=scriptname+'.Collection_make_sub: '

    log.info(logid+'STARTING PROCESSING FOR '+str(conditions))

    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    if combinations:
        combname = get_combo_name(combinations)
        for condition in combname:
            worklist = combname[condition]['works']
            envlist  = combname[condition]['envs']
            log.debug(logid+'LISTS: '+str(condition)+'\t'+str(worklist)+'\t'+str(envlist))
            subconf = NestedDefaultDict()
            add = list()

            smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'header.smk'))
            with open(smkf,'r') as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                    line = re.sub(condapath,'conda:  "../', line)
                    add.append(line)

            for i in range(len(worklist)):
                log.debug(logid+' SUBLISTS: '+str(worklist[i])+'\t'+str(envlist[i]))
                works = worklist[i].split('-')
                envs = envlist[i].split('-')
                subjobs = list()
                # Add variable for combination string
                subjobs.append('\ncombo = \''+str(envlist[i])+os.sep+'\'\n\nwildcard_constraints:\n\tcombo = combo,\n \trawfile = \'|\'.join(list(SAMPLES)),\n\tfile = \'|\'.join(list(samplecond(SAMPLES, config))),\n\tread = "R1|R2"\n')
                subjobs.append('\n\n')
                # Add rulethemall based on chosen workflows
                subjobs.append(''.join(rulethemall(subworkflows, config, loglevel, condapath, logfix, envlist[i])))        # RuleThemAll for snakemake depending on chosen workflows

                for j in range(len(works)):
                    listoftools, listofconfigs = create_subworkflow(config, works[j], [condition])

                    if listoftools is None:
                        log.warning(logid+'No entry fits condition '+str(condition)+' for processing step '+str(subwork))
                        return None

                    sconf = listofconfigs[0]
                    for a in range(0, len(listoftools)):
                        toolenv, toolbin = map(str, listoftools[a])
                        if toolenv != envs[j] or toolbin is None:
                            continue
                        sconf[works[j]+'ENV'] = toolenv
                        sconf[works[j]+'BIN'] = toolbin
                        subconf.update(sconf)
                        subname = toolenv+'.smk'

                        if works[j] == 'QC' and 'TRIMMING' in works and not 'MAPPING' in works:
                            if 'DEDUP' in works:
                                subname = toolenv+'_dedup_trim.smk'
                            else:
                                subname = toolenv+'_trim.smk'

                        if works[j] == 'QC' and not 'TRIMMING' in works and not 'MAPPING' in works:
                            if 'DEDUP' in subworkflows:
                                subname = toolenv+'_dedup.smk'
                            else:
                                subname = toolenv+'_raw.smk'

                        smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', subname))
                        with open(smkf,'r') as smk:
                            for line in smk.readlines():
                                line = re.sub(condapath, 'conda:  "../', line)
                                subjobs.append(line)
                            subjobs.append('\n\n')

                if 'MAPPING' in works:
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','mapping.smk'))
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')
                if 'QC' in subworkflows:
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','multiqc.smk'))
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

                # Append footer and write out subsnake and subconf per condition
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.smk'))
                with open(smkf,'r') as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath,'conda:  "../', line)
                        subjobs.append(line)
                    subjobs.append('\n\n')

                smko = os.path.abspath(os.path.join(subdir, '_'.join(['_'.join(condition), envlist[i], 'subsnake.smk'])))
                if os.path.exists(smko):
                    os.rename(smko, smko+'.bak')
                with open(smko, 'w') as smkout:
                    smkout.write(''.join(add))
                    smkout.write(''.join(subjobs))

                confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition), envlist[i], 'subconfig.json'])))
                if os.path.exists(confo):
                    os.rename(confo, confo+'.bak')
                with open(confo, 'w') as confout:
                    json.dump(subconf, confout)

                jobs.append([smko, confo])

    else:
        for condition in conditions:
            for subwork in subworkflows:
                log.debug(logid+'PREPARING '+str(subwork)+' '+str(condition))
                listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                if listoftools is None:
                    log.warning(logid+'No entry fits condition '+str(condition)+' for processing step '+str(subwork))
                    return None
                sconf = listofconfigs[0]
                for i in range(0, len(listoftools)):
                    toolenv, toolbin = map(str, listoftools[i])
                    if toolenv is None or toolbin is None:
                        continue
                    sconf[subwork+'ENV'] = toolenv
                    sconf[subwork+'BIN'] = toolbin
                    subconf.update(sconf)
                    subname = toolenv+'.smk'

                    add = ''.join(rulethemall(subworkflows))        # RuleThemAll for snakemake depending on chosen workflows

                    if subwork == 'QC' and 'TRIMMING' in workflow and not 'MAPPING' in workflow:
                        if 'DEDUP' in workflow:
                            subname = toolenv+'_dedup_trim.smk'
                        else:
                            subname = toolenv+'_trim.smk'

                    if subwork == 'QC' and not 'TRIMMING' in workflow and not 'MAPPING' in workflow:
                        if 'DEDUP' in subworkflows:
                            subname = toolenv+'_dedup.smk'
                        else:
                            subname = toolenv+'_raw.smk'

                    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', subname))
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

            if 'MAPPING' in works:
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','mapping.smk'))
                with open(smkf,'r') as smk:
                    for line in smk.readlines():
                        line = re.sub(condapath, 'conda:  "../', line)
                        subjobs.append(line)
                    subjobs.append('\n\n')
                if 'QC' in subworkflows:
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','multiqc.smk'))
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

            smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.smk'))
            with open(smkf,'r') as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                    line = re.sub(condapath,'conda:  "../', line)
                    subjobs.append(line)
                subjobs.append('\n\n')

            smko = os.path.abspath(os.path.join(tmpdir, '_'.join(['_'.join(condition), 'subsnake.smk'])))
            if os.path.exists(smko):
                os.rename(smko, smko+'.bak')
            with open(smko, 'a') as smkout:
                smkout.write(add)
                smkout.write(join('', subjobs))

            confo = os.path.abspath(os.path.join(tmpdir,'_'.join(['_'.join(condition), 'subconfig.json'])))
            if os.path.exists(confo):
                os.rename(confo, confo+'.bak')
            with open(confo, 'a') as confout:
                json.dump(subconf, confout)

            jobs.append([smko, confo])

    return jobs


@check_run
def make_post(postworkflow, config, samples, conditions, subdir, loglevel, subname=None, combinations=None):
    logid=scriptname+'.Collection_make_post: '

    log.info(logid+'STARTING POSTPROCESSING FOR '+str(conditions))

    jobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')
    summary_tools_set = set()
    summary_tools_dict = dict()

    if combinations:
        combname = get_combo_name(combinations)
        for condition in combname:
            subwork = postworkflow
            envlist = combname[condition]['envs']
            log.debug(logid+'POSTLISTS: '+str(condition)+'\t'+str(subwork)+'\t'+str(envlist))
            subconf = NestedDefaultDict()
            add = list()

            smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'header.smk'))
            with open(smkf, 'r') as smk:
                for line in smk.readlines():
                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                    line = re.sub(condapath,'conda:  "../', line)
                    add.append(line)

            for i in range(len(envlist)):
                envs = envlist[i].split('-')

                log.debug(logid+'POSTPROCESS: '+str(subwork)+' CONDITION: '+str(condition))
                if subwork in ['PEAKS', 'DE', 'DEU', 'DAS', 'DTU']:
                    listoftools, listofconfigs = create_subworkflow(config, subwork, combname)
                else:
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])

                if listoftools is None:
                    log.warning(logid+'No entry fits condition '+str(condition)+' for processing step '+str(subwork))
                    continue

                sconf = listofconfigs[0]
                for a in range(0, len(listoftools)):
                    subjobs = list()

                    toolenv, toolbin = map(str, listoftools[a])
                    if subwork in ['DE', 'DEU', 'DAS', 'DTU'] and toolbin not in ['deseq', 'diego']:  # for all other postprocessing tools we have more than one defined subworkflow
                        toolenv = toolenv+'_'+subwork

                    sconf[subwork+'ENV'] = toolenv
                    sconf[subwork+'BIN'] = toolbin

                    scombo = str(envlist[i])+os.sep if envlist[i] != '' else ''
                    combo = str.join(os.sep, [str(envlist[i]), toolenv])+os.sep if envlist[i] != '' else toolenv+os.sep

                    # Add variable for combination string
                    subjobs.append('\ncombo = \''+combo+'\'\n'+'\nscombo = \''+scombo+'\'\n'+'\nwildcard_constraints:\n\tcombo = combo,\n\tscombo = scombo,\n\tread = "R1|R2",\n\ttype = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"')
                    subjobs.append('\n\n')
                    subconf.update(sconf)

                    subname = toolenv+'.smk'
                    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', subname))

                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath, 'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

                    # Append footer and write out subsnake and subconf per condition
                    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'footer.smk'))
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(condapath,'conda:  "../', line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

                    smko = os.path.abspath(os.path.join(subdir, '_'.join(['_'.join(condition), envlist[i], subwork, toolenv, 'subsnake.smk'])))
                    if os.path.exists(smko):
                        os.rename(smko, smko+'.bak')
                    with open(smko, 'w') as smkout:
                        smkout.write(''.join(add))
                        smkout.write(''.join(subjobs))

                    confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition), envlist[i], subwork, toolenv, 'subconfig.json'])))
                    if os.path.exists(confo):
                        os.rename(confo, confo+'.bak')
                    with open(confo, 'w') as confout:
                        json.dump(subconf, confout)

                    jobs.append([smko, confo])

    return jobs


@check_run
def make_summary(summary_tools_set, summary_tools_dict, config, conditions, subdir, loglevel, subname=None, combinations=None):  # Need to check what we really need here, definitely conditions, config and workflows that should be summarized to retrieve combinations
    logid=scriptname+'.Collection_make_summary: '

    log.info(logid+'CREATING SUMMARY FOR '+str(conditions))

    jobs = list()
    subjobs = list()
    condapath = re.compile(r'conda:\s+"')
    logfix = re.compile(r'loglevel="INFO"')

    makeoutdir('REPORTS/SUMMARY')
    sum_path = os.path.join('nextsnakes', 'scripts', 'Analysis', 'SUMMARY')
    rmd_header = os.path.abspath(os.path.join(sum_path, 'header_summary.Rmd'))
    rmd_summary = os.path.abspath(os.path.join('REPORTS', 'SUMMARY', 'summary.Rmd'))
    combinations = get_combo(summary_tools_dict.keys(), config, conditions)

    if os.path.exists(rmd_summary):
        os.rename(rmd_summary, rmd_summary+'.bak')

    with open(rmd_header,'r') as read_file:
        for line in read_file.readlines():
            subjobs.append(line)
        subjobs.append('\n\n')

    for file in os.listdir(sum_path):
        file_path = os.path.abspath(os.path.join('nextsnakes', 'scripts', 'Analysis', 'SUMMARY', file))
        if file.startswith('SUM_'):
            file_tools = re.findall('[A-Z]+-[a-z]+', file)
            works = re.findall('[A-Z]+(?=-)', file)
            if set(file_tools).issubset(summary_tools_set):
                with open(file_path, 'r') as read_file:
                    if 'percompare' in os.path.basename(file_path):
                        for comparison in config[works[0]]['COMPARABLE'].keys():
                            for line in read_file.readlines():
                                line = re.sub('COMPARISON', comparison, line)
                                subjobs.append(line)
                            subjobs.append('\n\n')
                    else:
                        for line in read_file.readlines():
                            line = re.sub('COMPARISON', comparison, line)
                            subjobs.append(line)
                        subjobs.append('\n\n')

    with open(rmd_summary, 'a') as sumf:
        sumf.write(''.join(subjobs))
        sumf.write('\n\n')

    subjobs = list()

    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'header.smk'))
    with open(smkf,'r') as smk:
        for line in smk.readlines():
            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
            line = re.sub(condapath, 'conda:  "../', line)
            subjobs.append(line)
        subjobs.append('\n\n')

    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'summary.smk'))
    with open(smkf,'r') as smk:
        for line in smk.readlines():
            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
            line = re.sub(condapath, 'conda:  "../', line)
            subjobs.append(line)
        subjobs.append('\n\n')

    smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'footer.smk'))
    with open(smkf, 'r') as smk:
        for line in smk.readlines():
            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
            line = re.sub(condapath, 'conda:  "../', line)
            subjobs.append(line)
        subjobs.append('\n\n')

    smko = os.path.abspath(os.path.join(subdir, 'summary_subsnake.smk'))
    if os.path.exists(smko):
        os.rename(smko, smko+'.bak')
    with open(smko, 'a') as smkout:
        smkout.write(''.join(subjobs))
        smkout.write('\n\n')

    subconf = NestedDefaultDict()
    for key in ['BINS', 'MAXTHREADS', 'SETTINGS']:
        subconf[key] = config[key]
    subconf['WORKFLOWS'].merge(summary_tools_dict)

    confo = os.path.abspath(os.path.join(subdir, 'summary_subconfig.json'))
    if os.path.exists(confo):
        os.rename(confo, confo+'.bak')
    with open(confo, 'a') as confout:
        json.dump(subconf, confout)

    jobs.append([smko, confo])

    return jobs


@check_run
def rulethemall(subworkflows, config, loglevel, condapath, logfix, combo=''):
    logid=scriptname+'.Collection_rulethemall: '

    allmap = 'rule themall:\n\tinput:\texpand("MAPPED/{combo}{file}_mapped_sorted_unique.bam", combo=combo, file=samplecond(SAMPLES, config))' if not 'DEDUP' in subworkflows else 'rule themall:\n\tinput:\texpand("MAPPED/{combo}{file}_mapped_sorted_unique_dedup.bam", combo=combo, file=samplecond(SAMPLES, config))'
    allqc  = 'expand("QC/Multi/{combo}{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    allrawqc  = 'rule themall:\n\tinput:\texpand("QC/Multi/{combo}{condition}/multiqc_raw_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    alltrimqc = 'rule themall:\n\tinput:\texpand("QC/Multi/{combo}{condition}/multiqc_trimmed_raw_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'
    alltrim = 'rule themall:\n\tinput: expand("TRIMMED_FASTQ/{combo}{file}_{read}_trimmed.fastq.gz", combo=combo, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("TRIMMED_FASTQ/{combo}{file}_trimmed.fastq.gz", combo=combo, file=samplecond(SAMPLES, config))'
    alldedupqc = 'rule themall:\n\tinput:\texpand("QC/Multi/{combo}{condition}/multiqc_dedup_report.html", combo=combo, condition=str.join(os.sep, conditiononly(SAMPLES[0], config)))'
    alldedup = 'rule themall:\n\tinput: expand("DEDUP_FASTQ/{combo}{file}_{read}_dedup.fastq.gz", combo=combo, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("DEDUP_FASTQ/{combo}{file}_dedup.fastq.gz", combo=combo, file=samplecond(SAMPLES, config))'
    alltrimdedupqc = 'rule themall:\n\tinput:\texpand("QC/Multi/{combo}{condition}/multiqc_trim_dedup_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)'

    todos = list()

    if 'QC' in subworkflows and 'QC' in config:
        makeoutdir('QC')
        if 'MAPPING' in subworkflows:
            todos.append(allmap+',\n\t\t'+allqc+'\n\n')
        else:
            if 'TRIMMING' in subworkflows and 'DEDUP' not in subworkflows:
                todos.append(alltrimqc+'\n\n')
            elif 'TRIMMING' in subworkflows and 'DEDUP' in subworkflows:
                todos.append(alltrimdedupqc+'\n\n')
            elif 'DEDUP' in subworkflows and 'TRIMMING' not in subworkflows:
                todos.append(alldedupqc+'\n\n')
            else:
                todos.append(allrawqc+'\n\n')

    if 'MAPPING' in subworkflows and 'QC' not in subworkflows:
        log.info(logid+'Mapping without QC!')
        todos.append(allmap+'\n\n')

    if 'MAPPING' in subworkflows and 'TRIMMING' not in subworkflows:
        log.info(logid+'Simulated read trimming only!')
        makeoutdir('TRIMMED_FASTQ')
        smkf = os.path.abspath(os.path.join('nextsnakes', 'workflows', 'simulatetrim.smk'))
        with open(smkf,'r') as smk:
            todos.append(re.sub(condapath, 'conda:  "../', smk.read()))
        todos.append('\n\n')

    if 'TRIMMING' in subworkflows and 'QC' not in subworkflows and 'MAPPING' not in subworkflows:
        log.info(logid+'Trimming without QC!')
        todos.append(alltrim+'\n\n')

    if 'DEDUP' in subworkflows and 'QC' not in subworkflows and 'TRIMMING' not in subworkflows and 'MAPPING' not in subworkflows:
        log.info(logid+'DEDUP without QC!')
        todos.append(alldedup+'\n\n')

    return todos


@check_run
def tool_params(sample, runstate, config, subconf, tool = None):
    logid=scriptname+'.Collection_tool_params: '
    log.debug(logid+'Samples: '+str(sample))
    mp = OrderedDict()
    x = sample.split(os.sep)[:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    log.debug(logid+str([sample, runstate, subconf, x]))
    mp = subDict(config[subconf], x)[tool] if tool else subDict(config[subconf], x)
    log.debug(logid+'DONE: '+str(mp))
    return mp


@check_run
def get_reps(samples, config, analysis):
    logid=scriptname+'.Collection_get_reps: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        scond = sample.split(os.sep)[2:-1]
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(scond))
        partconf = subDict(config[analysis], scond)
        log.debug(logid+'CONF: '+str(partconf))
        ret['reps'].append(sample)
        wcfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        idx = partconf['REPLICATES'].index(wcfile)
        ret['pairs'].append(checkpaired_rep([str.join(os.sep, sample.split(os.sep)[2:])], config))
        ret['conds'].append(partconf['GROUPS'][idx])
        ret['batches'].append(partconf['BATCHES'][idx])
        if 'TYPES' in partconf and len(partconf['TYPES']) >= idx:
                ret['types'].append(str(partconf['TYPES'][idx]).replace(',','_'))
        else:
            ret['types'].append('std')

    rets = '-r '+str.join(',', ret['reps'])
    rets += ' -c '+str.join(',', ret['conds'])
    rets += ' -t '+str.join(',', ret['types'])
    rets += ' -b '+str.join(',', ret['batches'])
    rets += ' --paired '+str.join(',', ret['pairs']) if 'pairs' in ret else ''

    log.debug(logid+'RETURN: '+str(rets))
    return rets


@check_run
def get_diego_samples(samples, config, analysis):
    logid=scriptname+'.Collection_get_diego_samples: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(sample.split(os.sep)[1:-1]))
        partconf = subDict(config[analysis], sample.split(os.sep)[1:-1])
        log.debug(logid+'CONF: '+str(partconf))
        wcfile = str.join('-', sample.split(os.sep)[-4:]).replace('_mapped_sorted_unique.counts','')
        ret[wcfile].append(sample)

    log.debug(logid+'RETURN: '+str(ret))

    slist = ''
    for key, val in ret.items():
        slist +=  str(key)+'\t'
        slist += '\t'.join(val)
        slist += os.linesep

        log.debug(logid+'RETURN: '+str(slist))
    return slist


@check_run
def get_diego_groups(samples, config, analysis):
    logid=scriptname+'.Collection_get_diego_groups: '
    log.debug(logid+'Samples: '+str(samples))
    ret = defaultdict(list)
    for sample in samples:
        log.debug(logid+'WORKING ON: '+str(sample)+' CONDITION: '+str(sample.split(os.sep)[1:-1]))
        partconf = subDict(config[analysis], sample.split(os.sep)[1:-1])
        log.debug(logid+'CONF: '+str(partconf))
        wcfile = str.join('-', sample.split(os.sep)[-4:]).replace('_mapped_sorted_unique.counts','')
        checkfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        #wcfile = sample.split(os.sep)[-1].replace('_mapped_sorted_unique.counts','')
        idx = partconf['REPLICATES'].index(checkfile)
        cond = partconf['GROUPS'][idx]
        ret[cond].append(wcfile)

    slist = ''
    for key, val in ret.items():
        slist +=  str(key)+'\t'
        slist += '\t'.join(val)
        slist += os.linesep
    log.debug(logid+'RETURN: '+str(slist))
    return slist


@check_run
def env_bin_from_config(samples, config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config: '
    s = samples[0].split(os.sep)[:-1]
    mb, me = [None, None]
    for k in getFromDict(config[subconf], s):
        mb, me = k['BIN'], k['ENV']
    return mb, me


@check_run
def env_bin_from_config2(samples, config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config2: '
    for s in samples:
        log.debug(logid+'S: '+str(s))
        log.debug(logid+'C: '+str(conditiononly(s, config)))
        check = conditiononly(s, config)

        for k in getFromDict(config[subconf], check):
            if 'BIN' in k:
                mb = k['BIN']
            else:
                mb = ''
            if 'ENV' in k:
                me = k['ENV']
            else:
                me = ''
        log.debug(logid+str([str(mb), str(me)]))
    return mb, me


@check_run
def env_bin_from_config3(config, subconf):
    logid=scriptname+'.Collection_env_bin_from_config3: '
    envkey = subconf+"ENV"
    binkey = subconf+"BIN"
    me = config[envkey]
    mb = config[binkey]
    log.debug(logid+str([str(mb), str(me)]))
    return mb, me


@check_run
def rmempty(check):
    ret = list()
    for f in check:
        if os.path.isfile(f):
            ret.append(f)
    return ret


@check_run
def sample_from_path(path):
    ret = str(os.path.join(os.path.split(str(path))[-1]))
    return ret


@check_run
def runstate_from_sample(sample, config):
    logid = scriptname+'.Collection_runstate_from_sample: '
    ret = list()
    for s in sample:
        if len(getFromDict(config["SETTINGS"], s.split(os.sep))) < 1:
            s = os.path.dirname(s)
        n = s.split(os.sep)[-1]
        log.debug(logid+'SAMPLE: '+s)
        c = getFromDict(config["SETTINGS"], s.split(os.sep))[0]
        log.debug(logid+'SETTINGS: '+str(c))
        if dict_inst(c):
            if not c.get('SAMPLES'):
                for k, v in c.items():
                    log.debug(logid+'k, v: '+str([str(k), str(v)]))
                    #if n in v:
                    if k not in ret:
                        ret.append(k)
            else:
                ret.extend(s.split(os.sep))
        else:
            if n in c:
                k = s.split(os.sep)[-1]
                ret.append(k)
    log.debug(logid+'RETURN: '+str(ret))
    return ret


@check_run
def samplecond(sample, config): # takes list of sample names (including .fastq.gz) and returns a list with their conditions as directory path without fastq.gz ( ["condition/of/sample", ... ])
    logid = scriptname+'.Collection_samplecond: '
    ret = list()
    for s in sample:
        log.debug(logid+str(s))
        s = s.replace('.fastq.gz', '')
        check = os.path.dirname(s).split(os.sep)
        log.debug(logid+'CHECK: '+str(check))
        tmplist = check
        for r in runstate_from_sample([s], config):
            if r not in tmplist:
                tmplist.append(r)
        log.debug(logid+'TMPLIST: '+str(tmplist))
        paired = checkpaired([s], config)
        if 'paired' in paired:  # subDict(config['SETTINGS'], tmplist)['SEQUENCING']:
            s = re.sub(r'_[r|R|][1|2]', '', s)
        r = os.sep.join(tmplist)
        if r not in s:
            ret.append(os.sep.join([r, os.path.basename(s)]))
        else:
            ret.append(os.sep.join([os.path.dirname(s), os.path.basename(s)]))
    log.debug(logid+'RETURN: '+str(ret))
    return ret


@check_run
def conditiononly(sample, config):
    logid = scriptname+'.Collection_conditiononly: '
    ret = list()
    check = os.path.dirname(sample).split(os.sep)
    ret.extend(check)
    log.debug(logid+'CHECK: '+str(check))
    for r in runstate_from_sample([sample], config):
        log.debug(logid+'runstate '+str(r))
        if r not in ret:
            ret.append(r)
    log.debug(logid+'ret: '+str(ret))
    return ret


@check_run
def checkpaired(sample, config):
    logid = scriptname+'.Collection_checkpaired: '
    paired = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(s))
        check = conditiononly(s, config)
        log.debug(logid+'CHECK: '+str(check))
        p = subDict(config['SETTINGS'], check)
        paired = p.get('SEQUENCING')
        # Per sample paired, not implemented yet
        #pairedlist = p.get('SEQUENCING')
        #samplelist = p.get('SAMPLES')
        #x = samplelist.index(s.split(os.sep)[-1])
        #paired = pairedlist[x]
    log.debug(logid+'SEQUENCING: '+str(paired))
    return paired


@check_run
def checkpaired_rep(sample, config):
    logid = scriptname+'.Collection_checkpaired_rep: '
    log.debug(logid+'SAMPLE: '+str(sample))
    ret = list()
    for s in sample:
        check = conditiononly(s, config)
        p = subDict(config['SETTINGS'], check)
        paired = p.get('SEQUENCING')
        # Per sample paired, not implemented yet
        #pairedlist = p.get('SEQUENCING')
        #samplelist = p.get('SAMPLES')
        #x = samplelist.index(s.split(os.sep)[-1])
        #paired = pairedlist[x]
        ret.append(str(paired).replace(',','_'))
    log.debug(logid+'PAIRED: '+str(ret))
    return str.join(',', ret)


@check_run
def checkstranded(sample, config):
    logid = scriptname+'.Collection_checkstranded: '
    ret = list()
    stranded = ''
    for s in sample:
        check = conditiononly(s, config)
        p = subDict(config['SETTINGS'], check)
        log.debug(logid+'P: '+str(p))
        paired = p.get('SEQUENCING')
        # Per sample paired, not implemented yet
        #pairedlist = p.get('SEQUENCING')
        #samplelist = p.get('SAMPLES')
        #x = samplelist.index(s.split(os.sep)[-1])
        #paired = pairedlist[x]
        stranded = paired.split(',')[1] if len(paired.split(',')) > 1 else ''
    log.debug(logid+'STRANDEDNESS: '+str(stranded))
    return stranded


@check_run
def set_pairings(samples, config):
    logid = scriptname+'.Collection_set_pairings: '
    ret = list()
    log.debug(logid+'SAMPLES: '+str(samples))
    pairlist = config['PEAKS'].get('COMPARABLE')
    log.debug(logid+'PAIRLIST: '+str(pairlist))
    if pairlist:
        for k, v in pairlist.items():
            for x in samples:
                if k in x:
                    ret.extend(samplecond([x], config))
    else:
        return samples
    log.debug(logid+'return: '+str(ret))
    return ret


@check_run
def get_pairing(sample, stype, config, samples, scombo=''):
    logid = scriptname+'.Collection_get_pairings: '
    pairlist = config['PEAKS'].get('COMPARABLE')
    matching = ''
    log.debug(logid+'PAIRLIST: '+str(pairlist)+' SAMPLE: '+str(sample)+' SAMPLES: '+str(samples)+' Combo: '+str(scombo))
    if pairlist:
        for k, v in pairlist.items():
            if k in sample:
                for x in samples:
                    if v in x:
                        log.debug(logid+'Match found: '+str(v)+' : '+str(x))
                        matching = samplecond([x], config)[0].replace('MAPPED/', '')
                        log.debug(logid+'PAIRINGS: '+sample+': '+str(matching))
        log.debug(logid+'-c MAPPED/'+str(scombo)+str(matching)+'_mapped_'+str(stype)+'.bam')
        return '-c MAPPED/'+str(scombo)+str(matching)+'_mapped_'+str(stype)+'.bam'
    else:
        log.debug(logid+'No matching sample found')
        return ''


@check_run
def post_checkpaired(sample, config):
    logid = scriptname+'.Collection_checkpaired: '
    ret = list()
    paired = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(sample))
        check = conditiononly(sample, config)
        p = subDict(config['SETTINGS'], check)
        log.debug(logid+'P: '+str(p))
        paired = p.get('SEQUENCING').split(',')[0]
        #check = os.path.dirname(s).split(os.sep)
        #tmplist = check
        #p = getFromDict(config['SEQUENCING'], tmplist)[0]
        #if not dict_inst(p):
        #paired = p[0] if 'paired' in p or 'unpaired' in p or 'singlecell' in p else ''
        #for r in runstate_from_sample([s], config):
        #    log.debug(logid+'R: '+str(r))
        #    if r in p:
        #        tmplist.append(r)
        #        paired = getFromDict(config['SEQUENCING'], tmplist)[0].split(',')[0]
        #        tmplist = tmplist[:2]
    log.debug(logid+'PAIRED: '+str(paired))
    return paired


@check_run
def check_IP(sample, config):
    logid = scriptname+'.Collection_check_IP: '
    ret = list()
    clip = ''
    for s in sample:
        log.debug(logid+'SAMPLE: '+str(s))
        check = os.path.dirname(s).split(os.sep)
        r = runstate_from_sample([s], config)
        log.debug(logid+'RUNSTATE: '+str(r))
        tmplist = check
        if r not in tmplist:
            tmplist.extend(r)
        log.debug(logid+str(tmplist))
        log.debug(logid+'TMP: '+str(tmplist))
        check = getFromDict(config['PEAKS'], tmplist)[0]
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
        par = tool_params(sample, runstate , config, subconf)['OPTIONS'][idx]
        if par != '':
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
        log.debug(logid+'determine comparables in '+subwork)
        complist  = []
        compdict=config[subwork]['COMPARABLE']
        for key in compdict:
            for value in compdict[key]:
                complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr
    else:
        log.warning(logid+'no comparables found in '+subwork+'. Compare All vs. All.')
        groups_by_condition = list(yield_from_dict("GROUPS", config))
        flattened = sorted(set(val for sublist in groups_by_condition for val in sublist))
        combined = list(set(combinations(flattened,2)))
        complist = []
        for key, value in combined:
            complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr


@check_run
def comparable_as_string2(config, subwork):
    logid=scriptname+'.comparable_as_string2: '
    log.info(logid+'this is the config: '+str(config))
    check = config[subwork]['COMPARABLE']
    if check:
        log.debug(logid+'determine comparables in '+subwork)
        complist  = []
        compdict = config[subwork]['COMPARABLE']
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
        log.warning(logid+'no comparables found in '+subwork+'. Compare All vs. All.')
        groups_by_condition = list(yield_from_dict("GROUPS", config))
        flattened = sorted(set(val for sublist in groups_by_condition for val in sublist))
        combined = list(set(combinations(flattened,2)))
        complist = []
        for key, value in combined:
            complist.append(f"{key}vs{value}:{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr


@check_run
def comparable_as_string3(config, subwork):
    logid=scriptname+'.comparable_as_string: '
    check = config[subwork].get('COMPARABLE')
    if check:
        log.debug(logid+'determine comparables in '+subwork)
        complist  = []
        compdict=config[subwork]['COMPARABLE']
        for key in compdict:
            for value in compdict[key]:
                complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr
    else:
        log.warning(logid+'no comparables found in '+subwork+'. Compare All vs. All.')
        groups_by_condition = list(yield_from_dict("GROUPS", config))
        flattened = sorted(set(val for sublist in groups_by_condition for val in sublist))
        combined = list(set(combinations(flattened,2)))
        complist = []
        for key, value in combined:
            complist.append(f"{key}-vs-{value}")
        compstr = ','.join(complist)
        return compstr


##############################
########Nextflow Subs########
##############################
@check_run # SOURCE???
def nf_fetch_params(configfile):
    logid=scriptname+'.nf_fetch_params: '

    config = load_configfile(configfile)

    retconf = collections.defaultdict()
    retconf["REFERENCE"] = config["REFERENCE"]
    retconf["BINS"] = config["BINS"]
    retconf["MAXTHREAD"] = int(config["MAXTHREADS"])
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]
    retconf["SAMPLES"] = str.join(',', SAMPLES)
    LONGSAMPLES = samplecond(SAMPLES, config)
    retconf["LONGSAMPLES"] = str.join(',', LONGSAMPLES)
    log.info(logid+'Nextflow working on SAMPLES: '+str(SAMPLES))

    sample = SAMPLES[0]
    lsample = LONGSAMPLES[0]
    retconf["GENOME"] = genome(sample, config)
    retconf["SOURCE"] = source_from_sample(lsample, config)
    retconf["NAME"] = namefromfile(sample, config)
    paired = checkpaired([sample], config)
    retconf["PAIRED"] = paired
    stranded = checkstranded([sample], config)
    retconf["STRANDED"] = stranded

    if 'paired' in paired:
        log.info('RUNNING NEXTFLOW IN PAIRED READ MODE')
    if stranded != '':
        log.info('RUNNING NEXTFLOW WITH STRANDEDNESS '+str(stranded))

    if 'PEAKS' in config:
        retconf["CLIP"] = checkclip(SAMPLES, config)
        retconf["PEAKCONF"] = tool_params(sample, None, config,'PEAKS')['OPTIONS'][0]
        if 'ANNOTATION' in tool_params(sample, None, config,'PEAKS'):
            retconf["ANNOPEAK"] = tool_params(sample, None, config,'PEAKS')['ANNOTATION']
        else:
            retconf["ANNOPEAK"] = None
        try:
            all([x in peakconf for x in ['MINPEAKRATIO', 'PEAKDISTANCE', 'PEAKWIDTH', 'PEAKCUTOFF', 'MINPEAKHEIGHT', 'USRLIMIT']])
        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
            )
            log.error('Not all required peak finding options defined in config!\n'+''.join(tbe.format()))

        retconf["MINPEAKRATIO"] = peakconf['MINPEAKRATIO']
        retconf["PEAKDISTANCE"] = peakconf['PEAKDISTANCE']
        retconf["PEAKWIDTH"] = peakconf['PEAKWIDTH']
        retconf["PEAKCUTOFF"] = peakconf['PEAKCUTOFF']
        retconf["MINPEAKHEIGHT"] = peakconf['MINPEAKHEIGHT']
        retconf["USRLIMIT"] = peakconf['USRLIMIT']

        if 'PREPROCESS' in peakconf:
            retconf["PREPROCESS"] = ' '.join("{!s} {!s}".format(key, val) for (key, val) in peakconf['PREPROCESS'].items())

    return retconf


@check_run
def nf_tool_params(sample, runstate, config, subwork, toolenv, toolbin, workflows=None, condition=None):
    logid=scriptname+'.nf_tool_params: '
    log.debug(logid+'Samples: '+str(sample))
    t = genome(sample, config)
    mp = OrderedDict()
    x = sample.split(os.sep)[:-1]
    if runstate is None:
        runstate = runstate_from_sample([sample], config)[0]
    if runstate not in x:
        x.append(runstate)
    log.debug(logid+str([sample, runstate, config, t, x]))
    tp = list()

    if not workflows:
        mp = subDict(config[subwork], x)['OPTIONS']
        tp.append("--"+subwork+"ENV "+toolenv+" --"+subwork+"BIN "+toolbin+' ')
        if len(mp) > 0:
            for idx in range(len(mp)):
                tp.append("--"+toolenv+"_params_"+str(idx)+' \''+' '.join("{!s} {!s}".format(key, val) for (key, val) in mp[idx].items())+'\'')
    else:
        for subwork in workflows:
            sd = subDict(config[subwork], condition)
            mp = sd['OPTIONS']
            #listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
            #for i in range(0, len(listoftools)):
            toolenv, toolbin = map(str,[sd['ENV'], sd['BIN']])
            tp.append("--"+subwork+"ENV "+toolenv+" --"+subwork+"BIN "+toolbin+' ')

            if subwork == 'MAPPING':
                for idx in range(0,2):
                    tp.append("--"+toolenv+"_params_"+str(idx)+' \''+' '.join("{!s} {!s}".format(key, val) for (key, val) in mp[idx].items())+' \'')

                sfile = str.join(os.sep, x)+os.sep+sample.split(os.sep)[-1]
                ref = config["REFERENCE"]
                gdir = source_from_sample(sfile, config)
                mapper = toolenv
                name = namefromfile(sfile, config) or ''
                gen = genome(sfile, config)
                ext = check_tool_params(sfile, None , config, 'MAPPING', 2)
                uni = gen+name+'_'+ext
                anno = sd['ANNOTATION'] if 'ANNOTATION' in sd else None

                log.debug(logid+'FINAL PARAMS: '+' '.join([ref, gdir, mapper, name, gen, ext, uni, str(anno)]))

                index   = os.path.abspath(str.join(os.sep,[ref, gdir, mapper, ext, uni, mapper])+'.idx')
                reffa   = os.path.abspath(str.join(os.sep,[ref, gdir, gen+name])+'.fa.gz')
                genpath = os.path.abspath(str.join(os.sep,[ref, gdir, mapper, ext, uni]))
                anno    = os.path.abspath(str.join(os.sep,[ref, gdir, anno])) if anno else None

                apstr = "--"+subwork+"IDX "+index+" --"+subwork+"REF "+reffa+" --"+subwork+"GEN "+genpath
                if anno:
                    apstr += " --"+subwork+"ANNO "+anno+' '
                tp.append(apstr)

            else:
                if len(mp) > 0:
                    for idx in range(len(mp)):
                        tp.append("--"+toolenv+"_params_"+str(idx)+' \''+' '.join("{!s} {!s}".format(key, val) for (key, val) in mp[idx].items())+' \'')

    log.debug(logid+'DONE: '+str(tp))
    return ' '.join(tp)


##############################
#########Python Subs##########
##############################
@check_run
def dict_inst(d):
    logid=scriptname+'.Collection_dict_inst: '
    if isinstance(d, dict) or isinstance(d, OrderedDict) or isinstance(d, defaultdict) or isinstance(d, NestedDefaultDict):
        return True

@check_run
def getFromDict(dataDict, mapList):
    logid = scriptname+'.Collection_getFromDict: '
    log.debug(logid+'MAPLIST: '+str(mapList)+'\tDict: '+str(dataDict))
    ret = dataDict
    for k in mapList:
        log.debug(logid+'k: '+str(k))
        if k in dataDict:
            dataDict = dataDict[k]
            log.debug(logid+'subdict: '+str(dataDict))
        else:
            return list([])
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
        else:
            log.debug(logid+'No k in dict')
            return None
    return ret

@check_run
def subSetDict(dataDict, mapList):
    logid = scriptname+'.Collection_subDict: '
    log.debug(logid+str(mapList))
    parse = subDict(dataDict, mapList)
    ret = {}
    nested_set(ret, mapList, parse)
    log.debug(logid+str(ret))
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
def keysets_from_dict(dictionary, search=None, original=None):  # Only works for equal depth keysets, needs tweaking for other use cases
    logid = scriptname+'.Collection_keysets_from_dict: '

    keylist = list()
    if dict_inst(dictionary):
        for k, v in keys_from_dict(dictionary, search).items():
            keylist.append(v)
        log.debug(logid+'kl:'+str(keylist))
        combis = list()
        for i in range(1, len(keylist)+1):
            subkeylist = keylist[0:i]
            combis.extend(list(itertools.product(*subkeylist)))
        log.debug(logid+'cs:'+str(combis))
        ret = list()
        for combi in combis:
            #if len(getFromDict(dictionary, combi)) >= 1:
            check = subDict(dictionary, combi)
            log.debug(logid+'checking: '+str(check))
            if isvalid(check):
                if (isinstance(check, dict) and check.get("SAMPLES")) or isinstance(check, str):
                    log.debug(logid+'found: '+str(combi))
                    ret.append(combi)
        return ret
    else:
        return keylist

@check_run
def keys_from_dict(dictionary, search=None, save=None, first=True, lvl=0):
    logid = scriptname+'.Collection_keys_from_dict: '

    if first:
        first = False
        end = depth(dictionary)
        save = defaultdict(list)
        log.debug(logid+'TOTALDEPTH:'+str(end))

    if dict_inst(dictionary):
        log.debug(logid+'dictDEPTH: '+str(depth(dictionary)))
        log.debug(logid+'dictLEVEL: '+str(lvl))
        for k, v in dictionary.items():
            if not search or (search and k != search):
                save[lvl].append(k)
                log.debug(logid+'TMPSAVE: '+str(save))
                if dict_inst(v):
                    save = keys_from_dict(v, search, save, first, lvl+1)
                else:
                    continue
            else:
                log.debug(logid+'Found search: '+str(save))
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
        return dict(heapq.nsmallest(b, a.items(), key=itemgetter(1)))
    else:
        return dict({i:None for i in range(n)})

@check_run
def gethighest_dict(a, n):
    if n > len(a):
        b = len(a)
    else:
        b = n
    if len(a) > 0:
        return dict(heapq.nlargest(b, a.items(), key=itemgetter(1)))
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
        o.write(bytes(rand, encoding='UTF-8'))
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
        out += "%d\t%0.7f" % (it.index+1, it[0])+"\n"
        it.iternext()
    if o:
        o.write(bytes(out, encoding='UTF-8'))
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

@check_run
def multi_replace(repl, text):
    print('MULTI: '+str(repl)+str(text))
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, repl.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)

@check_run
def makelogdir(logdir):
    if not os.path.isabs(logdir):
        logdir =  os.path.abspath(logdir)
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    return logdir

@check_run
def get_dict_hash(d):
    logid = scriptname+'.get_dict_hash: '
    log.debug(logid+'INPUT DICT: '+str(d))
    ret = str(hashlib.sha256(bytes(str(sorted(d.items())),'utf-8')).hexdigest())
    log.debug(logid+'HASH: '+ret)
    return ret

########################################
############## DEPRECATED ##############
########################################

@check_run
def sources(config):
    logid = scriptname+'.Collection_sources: '
    ret = list()
    search =  [x[0] for x in keysets_from_dict(config["SOURCE"]) if x[0] != 'last']
    if len(getFromDict(config['SAMPLES'], search)) > 0:
        ret.extend(search)
    log.debug(logid+str(ret))
    return ret

@check_run
def genomepath(s, config):
    logid=scriptname+'.Collection_genomepath: '
    sa = os.path.basename(str(s))
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config["SAMPLES"])
    log.debug(logid+'GENOMEPATH: '+str([sa, cond, sk]))
    for skey in sk:
        klist = value_extract(skey, config["SOURCE"])
        for k in klist:
            for x, y in config["GENOME"].items():
                log.debug(logid+'GENOMEPATH: '+str([x, y, k]))
                if str(k) == str(y) or str(k) == str(x):
                    return os.path.join(str(x), str(y))

@check_run
def genome(s, config):
    logid=scriptname+'.Collection_genome: '
    sa = os.path.basename(str(s))
    sp = source_from_sample(str(s), config)
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
        ret.append(os.path.join(str(s), str(l)))
    return ret

@check_run
def genomename(s, config):
    s = os.path.basename(str(s))
    for k, v in config["SAMPLES"].items():
        for g, l in v.items():
            if s in l:
                for x, y in config["GENOME"].items():
                    if g == y:
                        return str(x)

@check_run
def transcriptomepath(s, config):
    logid=scriptname+'.Collection_transcriptomepath: '
    sa = os.path.basename(str(s))
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config["SAMPLES"])
    log.debug(logid+'TRANSCRIPTOMEPATH: '+str([sa, cond, sk]))
    for skey in sk:
        klist = value_extract(skey, config["SOURCE"])
        for k in klist:
            for x, y in config["TRANSCRIPTOME"].items():
                log.debug(logid+'TRANSCRIPTOMEPATH: '+str([x, y, k]))
                if str(k) == str(y) or str(k) == str(x):
                    return os.path.join(str(x), str(y))

@check_run
def transcriptome(s, config):
    logid=scriptname+'.Collection_transcriptome: '
    sa = os.path.basename(str(s))
    sp = source_from_sample(str(s), config)
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config['SAMPLES'])
    for skey in sk:
        klist = value_extract(skey, config['SOURCE'])
        for k in klist:
            for x, y in config['TRANSCRIPTOME'].items():
                log.debug(logid+str([k, x, y]))
                if str(k) == str(x):
                    return str(y)

@check_run
def fulltranscriptomepath(sa, config):
    ret=list()
    for s in sa:
        l = config["TRANSCRIPTOME"][s]
        ret.append(os.path.join(str(s), str(l)))
    return ret

@check_run
def transcriptomename(s, config):
    s = os.path.basename(str(s))
    for k, v in config["SAMPLES"].items():
        for g, l in v.items():
            if s in l:
                for x, y in config["TRANSCRIPTOME"].items():
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
def namefrompath(p, config):
    p = os.path.dirname(p).split(os.sep)
    klist = getFromDict(config["NAME"], p) if 'NAME' in config else list('')
    for k in klist:
        return str(k)

@check_run
def pathstogenomes(samples, config):
    ret = list()
    for s in samples:
        s = os.path.basename(s)
        for k, v in config["SAMPLES"].items():
            for g, l in v.items():
                if s in l:
                    for x, y in config["GENOME"].items():
                        if g == y:
                            ret.append(os.path.join(str(x), str(y)))
    return sorted(list(set(ret)))


@check_run
def source_from_sample(sample, config):
    logid=scriptname+'.Collection_source_from_sample: '
    s = os.path.dirname(str(sample))
    cond= s.split(os.sep)
    log.debug(logid+str([s, cond]))
    ret = getFromDict(config["SOURCE"], cond)[0]
    return ret


@check_run
def anno_from_file(sample, config, step):
    logid = scriptname+'.Collection_anno_from_file: '
    p = os.path.dirname(genomepath(sample, config))
    s = source_from_sample(sample, config)
    ret = os.path.join(config["REFERENCE"], p, subDict(config["ANNOTATE"], s)[step])
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
    ret = os.path.join(config["REFERENCE"], p, subDict(config["ANNOTATE"], lst)[step])
    log.debug(logid+str(ret))
    return ret

#
# Collection.py ends here
