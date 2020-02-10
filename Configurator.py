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
# Last-Updated: Mon Feb 10 17:16:45 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 293
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
    parser.add_argument("-f", "--preprocess", type=str, default='', help='Which preprocessing steps to conduct, choices are any or combinations of [\'SRA\', \'BASECALL\']. NOT IMPLEMENTED YET!!!')
    parser.add_argument("-w", "--workflows", type=str, default='', help='Which workflow steps to conduct, choices are any of or combinations of [\'MAPPING\', \'TRIMMING\', \'QC\']')
    parser.add_argument("-l", "--postprocess", type=str, default='', help='Which workflow steps to conduct,choices are any of or combinations of [\'COUNTING\',\'UCSC\',\'PEAKS\',\'ANNOTATE\',\'DE\',\'DEU\']')
    parser.add_argument("-g", "--gcs", type=str, default='homo_sapiens:condition:setting', help='Comma separated list of colon separated GenomeConditionSetting relationship. For each genome to work on you can define one or multiple conditions and settings that will be used for the analysis, e.g. hg38:WT:singleend,hg38:KO:pairedend,dm6:01012020:testsequencing or just a single colon separated GCS')
    parser.add_argument("-r", "--refdir", type=str, default='GENOMES', help='Path to directory with reference genome')
    parser.add_argument("-m", "--genomemap", type=str, default='homo_sapiens:hg38', help='Comma separated list of colon separated mapping of genome(s) of interest to genome FASTA.gz file, e.g. homo_sapiens:hg38,drosophila:dm6')
    parser.add_argument("-x", "--genomeext", type=str, default='homo_sapiens:None', help='Comma separated list of colon separated mapping of genome(s) of interest to extensio in FASTA.gz file, e.g. homo_sapiens:_extended,drosophila:minimal. This is not required if there is no extension which is often the case.')
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
def create_json_config(configfile, append, skeleton, preprocess, workflows, postprocess, gcs, refdir, binaries, procs, scripts, genomemap, genomeext, optionalargs=None):

    # CLEANUP
    oldcnf = os.path.abspath(configfile)
    for oldfile in glob.glob(oldcnf):
        shutil.copy2(oldfile,oldfile+'.bak')
        log.warning(logid+'Found old config file'+oldfile+' created backup of old config '+oldfile+'.bak')

    config = load_configfile(os.path.abspath(skeleton))
    newconf = NestedDefaultDict()

    todos = ','.join([x for x in [preprocess,workflows,postprocess] if x is not '' ]).split(',')
    for x in todos:
        if x not in config:
            log.error(logid+'Key '+str(x)+' not found in skeleton, please check for typos!')
            sys.exit()

    log.info(logid+'Creating config json for steps '+str(todos))

    if genomemap:
        genmap = [x.split(':') for x in genomemap.split(',')]
    else:
        if not append:
            log.error(logid+'No mapping of genome to genome fasta found, please provide -m option')
            sys.exit()

    if genomeext:
        genext = [x.split(':') for x in genomeext.split(',')]

    if not append:
        #newconf.merge(config)
        newconf['PREPROCESSING'] = preprocess
        newconf['WORKFLOWS'] = workflows
        newconf['POSTPROCESSING'] = postprocess
        newconf['REFERENCE'] = refdir
        newconf['BINS'] = binaries
        newconf['MAXTHREADS'] = str(procs)
        newconf['GENOME'] = NestedDefaultDict()
        for g in genmap:
            newconf['GENOME'][str(g[0])] = str(g[1])

    else:
        oldconfig = load_configfile(os.path.abspath(os.path.join(configfile)))
        #newconf.merge(oldconfig)

        if preprocess and preprocess not in newconf['PREPROCESSING']:
            newconf['PREPROCESSING'] = str.join(',',[oldconf['PREPROCESSING'],preprocess])
        if workflows and workflows not in newconf['WORKFLOWS']:
            newconf['WORKFLOWS'] = str.join(',',[oldconf['WORKFLOWS'],workflows])
        if postprocess and postprocess not in newconf['POSTPROCESSING']:
            newconf['POSTPROCESSING'] = str.join(',',[oldconf['POSTPROCESSING'],postprocess])
        if refdir and refdir != oldconf['REFERENCE']:
            newconf['REFERENCE'] = refdir
        else:
            newconf['REFERENCE'].merge(oldconf['REFERENCE'])
        if binaries and binaries != oldconf['BINS']:
            newconf['BINS'] = binaries
        else:
            newconf['BINS'].merge(oldconf['BINS'])
        if procs and procs != oldconf['MAXTHREADS']:
            newconf['MAXTHREADS'] = str(procs)
        else:
            newconf['MAXTHREADS'].merge(oldconf['MAXTHREADS'])
        if genomemap and genmap[0] not in newconf['GENOME'] or newconf['GENOME'][genmap[0]] != genmap[1]:
            newconf['GENOME'] = NestedDefaultDict()
            for g in genmap:
                newconf['GENOME'][str(g[0])] = str(g[1])
        else:
            newconf['GENOME'] = oldconf['GENOME']

        for key in ['NAME','SOURCE','SAMPLES','SEQUENCING']:
            for genome,condition,setting in [x.split(':') for x in gcslist]:
                if key == 'NAME':
                    if genomeext:
                        for x in genext:
                            if str(genome) == str(x[0]):
                                if str(x[1]) != oldconfig[key][genome][condition][setting]:
                                    newconf[key][genome][condition][setting] = str(x[1])
                                else:
                                    newconf[key][genome][condition][setting] = oldconf[key][genome][condition][setting]
                elif key == 'SOURCE':
                    if genomemap:
                        for x in genmap:
                            if str(genome) == str(x[0]) :
                                if str(x[1]) != oldconfig[key][genome][condition][setting]:
                                    newconf[key][genome][condition][setting] = str(x[1])
                                else:
                                    newconf[key][genome][condition][setting] = oldconf[key][genome][condition][setting]
                else:
                    newconf[key][genome][condition][setting] = oldconf[key][genome][condition][setting]

        for do in todos:
            if do not in newconf:
                newconf[do].merge(oldconf[do])

    """Now we replace the placeholders in the skeleton config with the actual ones or update an existing config with new workflows"""

    if gcs or append:
        if append:
            gcstemplist = list_all_keys_of_dict(oldconf['SAMPLES'])
            log.debug(logid+'OLDCONF: '+str(oldconf["SAMPLES"]))
            for g in gcstemplist:
                print(g)
            sys.exit()
        else:
            gcslist = gcs.split(',')
    else:
        log.error(logid+'GenomeConditionSetting (gcs) not defined!')
        sys.exit()

    log.debug(logid+'List of GenomeConditionSettings: '+str(gcs)+'\t'+str(gcslist))

    for key in ['NAME','SOURCE','SAMPLES','SEQUENCING']:
        for genome,condition,setting in [x.split(':') for x in gcslist]:
            if key == 'NAME':
                if genomeext and not append:
                    for x in genext:
                        if str(genome) == str(x[0]) :
                            if str(x[1]) is None or str(x[1]) == 'None':
                                x[1] = ''
                            newconf[key][genome][condition][setting] = str(x[1])
                else:
                    newconf[key][genome][condition][setting] = config[key]['genome']['condition']['setting']
            elif key == 'SOURCE':
                if genomemap:
                    for x in genmap:
                        if str(genome) == str(x[0]) :
                            newconf[key][genome][condition][setting] = str(x[1])
                else:
                    newconf[key][genome][condition][setting] = config[key]['genome']['condition']['setting']
            else:
                newconf[key][genome][condition][setting] = config[key]['genome']['condition']['setting']

    for do in todos:
        if do not in newconf:
            newconf[do].merge(config[do])

    for key in todos:
        log.debug(logid+'OLD: '+str(key)+'\t'+str(config[key]))
        for genome,condition,setting in [x.split(':') for x in gcslist]:
            if genome not in newconf[key]:
                newconf[key][genome] = NestedDefaultDict()
                log.debug(logid+'Genome: '+str(newconf[key]))
            if condition not in newconf[key][genome]:
                newconf[key][genome][condition] = NestedDefaultDict()
                log.debug(logid+'Condition: '+str(newconf[key]))
            if setting not in newconf[key][genome][condition]:
                newconf[key][genome][condition][setting] = NestedDefaultDict()
                log.debug(logid+'SETTING: '+str(newconf[key]))

            if 'genome' in newconf[key]:
                newconf[key][genome] = newconf[key].pop('genome')
                newconf[key][genome][condition] = newconf[key][genome].pop('condition')
                newconf[key][genome][condition][setting] = newconf[key][genome][condition].pop('setting')
            else:
                log.debug(logid+'TODO: '+str(key)+'\t'+str(config[key])+'\t'+str(newconf[key]))
                newconf[key][genome][condition][setting].merge(config[key]['genome']['condition']['setting'])

    print_json(newconf,configfile)

@check_run
def print_json(paramdict,ofn):
    with open(ofn,'w') as jsonout:
        json.dump(paramdict,jsonout)

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


        create_json_config(knownargs.configfile, knownargs.append, knownargs.skeleton, knownargs.preprocess, knownargs.workflows, knownargs.postprocess, knownargs.gcs, knownargs.refdir, knownargs.binaries, knownargs.procs, knownargs.scripts, knownargs.genomemap, knownargs.genomeext, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))



# Configurator.py ends here
