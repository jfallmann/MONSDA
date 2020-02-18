#!/usr/bin/env python3
# Configurator.py ---
#
# Filename: Configurator.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Feb 10 08:09:48 2020 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Feb 18 11:37:09 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 505
# URL:
# Doc URL:
# Keywords:
# Compatibility:
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
# <http://www.gnu.org/licenses/>.

import glob, os, sys, inspect, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re
min_version("5.8.2")

from lib.Collection import *
from lib.Logger import *
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Helper to create initial config file used for workflow processing')
    parser.add_argument("-c", "--configfile", type=str, default='configurator.json', help='Configuration json to write to, can be called together with --append option to append new workflows to existing config')
    parser.add_argument("-a", "--append", action="store_true", help='If set configuration will be appended to existing json')
    parser.add_argument("-s", "--skeleton", type=str, default='snakes/configs/skeleton.json', help='Skeleton config to build from, per default the one that comes with this repository, change only when you know what you do')
    parser.add_argument("-p", "--preprocess", type=str, default='', help='Which preprocessing steps to conduct, choices are any or combinations of [\'SRA\', \'BASECALL\']. NOT IMPLEMENTED YET!!!')
    parser.add_argument("-w", "--workflows", type=str, default='', help='Which workflow steps to conduct, choices are any of or combinations of [\'MAPPING\', \'TRIMMING\', \'QC\']')
    parser.add_argument("-l", "--postprocess", type=str, default='', help='Which workflow steps to conduct,choices are any of or combinations of [\'COUNTING\',\'UCSC\',\'PEAKS\',\'ANNOTATE\',\'DE\',\'DEU\']')
    parser.add_argument("-r", "--refdir", type=str, default='GENOMES', help='Path to directory with reference genome')
    parser.add_argument("-i", "--ics", type=str, default='id:condition:setting', help='Comma separated list of colon separated IdentifierConditionSetting relationship. For each id to work on you can define one or multiple conditions and settings that will be used for the analysis, e.g. hg38:WT:singleend,01012020:KO:pairedend,X321F5:01012020:testsequencing or just a single colon separated ICS')
    parser.add_argument("-g", "--genomes", type=str, default='hg38:hg38', help='Comma separated list of colon separated mapping of genome-IDs to genome FASTA.gz filename, e.g. hg38:hg38,01012020:dm6 means ID hg38 maps to a file hg38.fa.gz and ID 01012020 maps to a file dm6.fa.gz')
    parser.add_argument("-m", "--genomemap", type=str, default='id:hg38', help='Comma separated list of colon separated mapping of sample-IDs to genome-IDs, e.g. hg38:hg38,01012020:dm6')
    parser.add_argument("-x", "--genomeext", type=str, default=None, help='Comma separated list of colon separated mapping of genome-IDs to extension in FASTA.gz file, e.g. hg38:_extended,01012020:_bisulfit. This is not required if there is no extension which is often the case.')
    parser.add_argument("--binaries", type=str, default='snakes/scripts', help='Path to binary directory')
    parser.add_argument("-b", "--scripts", type=str, default='snakes/scripts', help='Path to script for execution')
    parser.add_argument("-j", "--procs", type=int, default=1, help='Maximum number of parallel processes to start snakemake with, represented by MAXTHREADS in config')
    parser.add_argument("-v", "--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

####################
#### FUNCTIONS  ####
####################

def check_run(func):
    def func_wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)

        except Exception as err:
            exc_type, exc_value, exc_tb = sys.exc_info()
            tbe = tb.TracebackException(
                exc_type, exc_value, exc_tb,
            )
            log.error(''.join(tbe.format()))
    return func_wrapper

@check_run
def create_json_config(configfile, append, skeleton, preprocess, workflows, postprocess, ics, refdir, binaries, procs, scripts, genomemap, genomes, genomeext, optionalargs=None):

    # CLEANUP
    oldcnf = os.path.abspath(configfile)
    for oldfile in glob.glob(oldcnf):
        shutil.copy2(oldfile,oldfile+'.bak')
        log.warning(logid+'Found old config file'+oldfile+' created backup of old config '+oldfile+'.bak')

    config = load_configfile(os.path.abspath(skeleton))
    newconf = NestedDefaultDict()
    oldconf = NestedDefaultDict()
    icslist = list()

    todos = ','.join([x for x in [preprocess,workflows,postprocess] if x is not '' ]).split(',')
    for x in todos:
        if x not in config:
            log.error(logid+'Key '+str(x)+' not found in skeleton, please check for typos!')
            sys.exit()

    log.info(logid+'Creating config json for steps '+str(todos))

    genmap = defaultdict()
    if genomemap:
        genmap = {key: value for (key, value) in [x.split(':') for x in genomemap.split(',')]}
        log.debug(logid+'GENOMEMAP: '+str(genmap))
    else:
        if not append:
            log.error(logid+'No mapping of sample-ID to genome-ID found, please provide -m option')
            sys.exit()

    gens = defaultdict()
    if genomes:
        gens = {key: value for (key, value) in [x.split(':') for x in genomes.split(',')]}
        log.debug(logid+'GENOMES: '+str(gens))
    else:
        if not append:
            log.error(logid+'No mapping of genome to genome fasta found, please provide -g option')
            sys.exit()

    genext = defaultdict()
    if genomeext:
        genext = {key: value for (key, value) in [x.split(':') for x in genomeext.split(',')]}
        log.debug(logid+'GENOMEEXTENSION: '+str(genext))
    if ics or append:
        if append:
            oldconf = load_configfile(os.path.abspath(os.path.join(configfile)))
            iteration = -1
            icstemp = ''
            for k,v in list_all_keys_of_dict(oldconf['SAMPLES']):
                iteration+=1
                if k == 'last':
                    icslist.append(icstemp[:-1])
                    if iteration >3:
                        icstemp = icstemp.split(':')[0]+':'
                        iteration = -1
                    else:
                        icstemp=''
                        iteration = -1
                else:
                    icstemp+=k+':'
            if ics:
                for x in ics.split(','):
                    if x not in icslist:
                        icslist.append(x)
        else:
            icslist = ics.split(',')
    else:
        log.error(logid+'IdentifierConditionSetting (ics) not defined!')
        sys.exit()

    log.debug(logid+'List of IdentifierConditionSettings: '+str(icslist))


    if not append:
        #newconf.merge(config)
        newconf['PREPROCESSING'] = preprocess
        newconf['WORKFLOWS'] = workflows
        newconf['POSTPROCESSING'] = postprocess
        newconf['REFERENCE'] = refdir
        newconf['BINS'] = binaries
        newconf['MAXTHREADS'] = str(procs)
        newconf['GENOME'] = NestedDefaultDict()
        for k,v in gens.items():
            newconf['GENOME'][str(k)] = str(v)

        for key in ['NAME','SOURCE','SEQUENCING','SAMPLES']:
            for id,condition,setting in [x.split(':') for x in icslist]:
                if key == 'NAME':
                    if genomeext:
                        for k,v in genext.items():
                            if str(v) is None or str(v) == 'None':
                                v = ''
                            newconf[key][id][condition][setting] = str(v)
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                elif key == 'SOURCE':
                    if genomemap:
                        for k,v in genmap.items():
                            if k == id:
                                newconf[key][id][condition][setting] = str(v)
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                elif key == 'SAMPLES':
                    samplelist = get_samples_from_dir(id, condition, setting, newconf)
                    if samplelist:
                        newconf[key][id][condition][setting] = samplelist
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                else:
                    newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']


    else:
        #newconf.merge(oldconfig)

        if preprocess and preprocess not in newconf['PREPROCESSING']:
            newconf['PREPROCESSING'] = str.join(',',list(set(str.join(',',[oldconf['PREPROCESSING'],preprocess]).split(','))))
        if workflows and workflows not in newconf['WORKFLOWS']:
            newconf['WORKFLOWS'] = str.join(',',list(set(str.join(',',[oldconf['WORKFLOWS'],workflows]).split(','))))
        if postprocess and postprocess not in newconf['POSTPROCESSING']:
            newconf['POSTPROCESSING'] = str.join(',',list(set(str.join(',',[oldconf['POSTPROCESSING'],postprocess]).split(','))))
        if refdir and refdir != oldconf['REFERENCE']:
            newconf['REFERENCE'] = refdir
        else:
            newconf['REFERENCE'] = str(oldconf['REFERENCE'])
        if binaries and binaries != oldconf['BINS']:
            newconf['BINS'] = binaries
        else:
            newconf['BINS'] = str(oldconf['BINS'])
        if procs and procs != oldconf['MAXTHREADS']:
            newconf['MAXTHREADS'] = str(procs)
        else:
            newconf['MAXTHREADS'] = str(oldconf['MAXTHREADS'])

        log.debug(logid+'GENOMEMAP: '+str(genomemap)+'\t'+str(genmap))
        if genomes and any([x not in newconf['GENOME'] for x in list(gens.keys())]) or any([[x not in newconf['GENOME'][y] for x in gens[y]] for y in gens.keys()]):
            newconf['GENOME'] = NestedDefaultDict()
            newconf['GENOME'].merge(oldconf['GENOME'])
            for k,v in gens.items():
                newconf['GENOME'][str(k)] = str(v)
        else:
            newconf['GENOME'] = str(oldconf['GENOME'])

        log.debug(logid+'GENOMEMAPCONF: '+str(newconf['GENOME']))

        for key in ['NAME','SOURCE','SAMPLES','SEQUENCING']:
            for id,condition,setting in [x.split(':') for x in icslist]:
                if key == 'NAME' or key == 'SOURCE':
                    try:
                        checkkey=getFromDict(oldconf[key],[id,condition,setting])
                    except:
                        checkkey=list()
                    if len(checkkey) > 0:
                        if key == 'NAME':
                            if genomeext:
                                for k,v in genext.items():
                                    if id in [x for x in find_key_for_value(k,genmap)]:
                                        if str(v) != oldconf[key][id][condition][setting]:
                                            newconf[key][id][condition][setting] = str(v)
                                        else:
                                            newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                            else:
                                newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                        elif key == 'SOURCE':
                            if genomemap:
                                for k,v in genmap.items():
                                    if k == id:
                                        if str(v) != oldconf[key][id][condition][setting]:
                                            newconf[key][id][condition][setting] = str(v)
                                        else:
                                            newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                            else:
                                newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                        else:
                            newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                else:
                    newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']

        for do in todos:
            if do not in newconf and do in oldconf:
                newconf[do].merge(oldconf[do])

    """Now we replace the placeholders in the skeleton config with the actual ones or update an existing config with new workflows"""

    log.debug(logid+'NEW: '+str(newconf))

    for do in todos:
        if do not in newconf:
            newconf[do].merge(config[do])

    for key in todos:
        log.debug(logid+'OLD: '+str(key)+'\t'+str(config[key]))
        for id,condition,setting in [x.split(':') for x in icslist]:
            if id not in newconf[key]:
                newconf[key][id] = NestedDefaultDict()
                log.debug(logid+'ID: '+str(newconf[key]))
            if condition not in newconf[key][id]:
                newconf[key][id][condition] = NestedDefaultDict()
                log.debug(logid+'Condition: '+str(newconf[key]))
            if setting not in newconf[key][id][condition]:
                newconf[key][id][condition][setting] = NestedDefaultDict()
                log.debug(logid+'SETTING: '+str(newconf[key]))

            if 'id' in newconf[key]:
                newconf[key][id] = newconf[key].pop('id')
                newconf[key][id][condition] = newconf[key][id].pop('condition')
                newconf[key][id][condition][setting] = newconf[key][id][condition].pop('setting')
            else:
                log.debug(logid+'TODO: '+str(key)+'\t'+str(config[key])+'\t'+str(newconf[key]))
                newconf[key][id][condition][setting].update(config[key]['id']['condition']['setting'])

    print_json(newconf,configfile)

@check_run
def print_json(paramdict,ofn):
    with open(ofn,'w') as jsonout:
        print(json.dumps(paramdict,indent=4),file=jsonout)

####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        knownargs=args[0]
        optionalargs=args[1:]

        makelogdir('LOGS')
        log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=knownargs.loglevel)
        log.addHandler(logging.StreamHandler(sys.stdout))  # streamlog

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')


        create_json_config(knownargs.configfile, knownargs.append, knownargs.skeleton, knownargs.preprocess, knownargs.workflows, knownargs.postprocess, knownargs.ics, knownargs.refdir, knownargs.binaries, knownargs.procs, knownargs.scripts, knownargs.genomemap, knownargs.genomes, knownargs.genomeext, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))



# Configurator.py ends here
