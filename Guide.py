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
# Last-Updated: Fri May  8 12:10:43 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 552
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
from collections import defaultdict, deque
import traceback as tb
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re
import sys
import copy
import json
#import logging
min_version("5.8.2")

from lib.Collection import *
from lib.Logger import *
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Helper to create initial config file used for workflow processing')
    parser.add_argument("-m", "--manualmode", action="store_true", help='If set configuration will explain every step in detail')

    return parser.parse_args()

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
def create_json_config(configfile, append, template, preprocess, workflows, postprocess, ics, refdir, binaries, procs, genomemap, genomes, genomeext, sequencing, annotation, optionalargs=None):
    # CLEANUP
    oldcnf = os.path.abspath(configfile)
    for oldfile in glob.glob(oldcnf):
        shutil.copy2(oldfile,oldfile+'.bak')
        log.warning(logid+'Found old config file'+oldfile+' created backup of old config '+oldfile+'.bak')

    config = load_configfile(os.path.abspath(template))
    newconf = NestedDefaultDict()
    oldconf = NestedDefaultDict()
    icslist = list()

    todos = ','.join([x for x in [preprocess,workflows,postprocess] if x != '' ]).split(',')
    for x in todos:
        if x not in config:
            log.error(logid+'Key '+str(x)+' not found in template, please check for typos!')
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

    seqlist = [s.replace(':',',') for s in sequencing.split(',')]

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
                            # if v in newconf['GENOME']:
                            newconf[key][id][condition][setting] = str(v)
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                elif key == 'SEQUENCING':
                    if len(seqlist) > 0:
                        newconf[key][id][condition][setting] = deque(seqlist).popleft()
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                elif key == 'SAMPLES':
                    samplelist = get_samples_from_dir(id, condition, setting, newconf)
                    log.debug(logid+'SAMPLELIST: '+str(samplelist))
                    if len(samplelist) > 0:
                        newconf[key][id][condition][setting] = samplelist
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                else:
                    newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']

    else:
        # newconf.merge(oldconfig)

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
                                    if v in newconf['GENOME']:
                                        if str(v) != str(oldconf[key][id][condition][setting]):
                                            newconf[key][id][condition][setting] = str(v)
                                        else:
                                            newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                            else:
                                newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                        else:
                            newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                elif key == 'SAMPLES':
                    samplelist = get_samples_from_dir(id, condition, setting, oldconf)
                    log.debug(logid+'SAMPLELIST: '+str(samplelist))
                    if len(samplelist) > 0:
                        newconf[key][id][condition][setting] = samplelist
                    else:
                        newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']
                else:
                    newconf[key][id][condition][setting] = config[key]['id']['condition']['setting']

        for do in todos:
            if do not in newconf and do in oldconf:
                newconf[do].merge(oldconf[do])

    """Now we replace the placeholders in the template config with the actual ones or update an existing config with new workflows"""

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

    # print_json(newconf,configfile,annotation)
    # input=json.dumps(newconf)
    # flatconf=json.loads(input)
    return newconf


# @check_run
def print_json(paramdict,ofn,annotation=''):
    with open(ofn,'w') as jsonout:
        if annotation != '' and 'gff.gz' in annotation or 'gff3.gz' in annotation:
            print(re.sub('genome_or_other.gff3.gz',annotation,json.dumps(paramdict,indent=4)),file=jsonout)
        else:
            print(json.dumps(paramdict,indent=4),file=jsonout)

def conversation(question, origin):
    print(question)
    default = copy.deepcopy(origin)
    if default is not None:
        print("\n")

        if isinstance(default, str):
            print("\t> ",default)
            print("\n")
            print("enter what should be added or enter 'n' to continue ")
            a = input(">>> ")
            if a == 'n':
                print("fine, everything the same!")
                return default
            else:
                return a

        if isinstance(default, dict):
            if default:
                for element in default:
                    print("\t> ",element,": ",default[element])
            else:
                print("\t> no entrys so far")
            print("\n")
            while True:
                print("enter 'y' for yes or 'n' for no")
                a = input(">>> ")
                if a=='n':
                    return
                if a=='y':
                    print("okay, tell me your setting:")
                    if default:
                        for element in default:
                            print(element+":")
                            a = input(">>> ").strip().replace(" ","")
                            default[element]=a
                        return default
                    else:
                        while True:
                            print("enter the name of the argument or enter 'n' to quit")
                            name = input(">>> ").strip().replace(" ","")
                            if name=='n':
                                break
                            print("enter the value")
                            value = input(">>> ").strip().replace(" ","")
                            default.update({name:value})

                    return default

    a = input(">>> ").strip().replace(" ","")
    return a

def rename(config,job,*args):
    ics_list,naming_list = ics_naming(config,job)
    for naming in naming_list:
        switch=False
        for ics in ics_list:
                if naming==ics:
                    switch=True
                    print(f"okay, look at the settings of {ics[0][0]}:{ics[1][0]}:{ics[2][0]}")
                    for call in args:
                        if len(call)==1:
                            text  =call[0]
                            config[job][ics[0][0]][ics[1][0]][ics[2][0]] = conversation(text,None)
                        if len(call)==2:
                            text  =call[0]
                            key   =call[1]
                            config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = conversation(text,None)
                        if len(call)==3:
                            text  =call[0]
                            key   =call[1]
                            number=call[2]
                            config[job][ics[0][0]][ics[1][0]][ics[2][0]][key][number] = conversation(text, config[job][naming[0][0]][naming[1][0]][naming[2][0]][key][number])
                elif ics in naming_list and naming!=ics:
                    switch=False
                else:
                    if switch:
                        for call in args:
                            if len(call)==1:
                                text  =call[0]
                                config[job][ics[0][0]][ics[1][0]][ics[2][0]] = config[job][naming[0][0]][naming[1][0]][naming[2][0]]
                            if len(call)==2:
                                text  =call[0]
                                key   =call[1]
                                config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = config[job][naming[0][0]][naming[1][0]][naming[2][0]][key]
                            if len(call)==3:
                                text  =call[0]
                                key   =call[1]
                                number=call[2]
                                config[job][ics[0][0]][ics[1][0]][ics[2][0]][key][number] = config[job][naming[0][0]][naming[1][0]][naming[2][0]][key][number]

def explain(text):
    if manual:
        print("\n** ",text,"\n")


def get_path(input_dict):
    for key, value in input_dict.items():
        if isinstance(value, dict):
            for subkey in get_path(value):
                yield key + ':' + subkey
        else:
            yield key

def ics_naming(config,workflow):
    ics_list=[]
    naming_list=[]
    level = conversation(f"set the following {workflow} settings per 'setting', per 'condition', per 'id', or 'all' the same?",None)
    a=0
    for id in config[workflow].keys():
        if id == "TOOLS" or id == "COMPARABLE":
            continue
        i=0
        for condition in config[workflow][id].keys():
            c=0
            for setting in config[workflow][id][condition].keys():
                ics_list.append([[id],[condition],[setting]])
                if level=='setting':
                    naming_list.append([[id],[condition],[setting]])
                if level=='condition' and c==0:
                    naming_list.append([[id],[condition],[setting]])
                    c+=1
                if level=='id' and i==0:
                    naming_list.append([[id],[condition],[setting]])
                    i+=1
                if level=='all' and a==0:
                    naming_list.append([[id],[condition],[setting]])
                    a+=1
    return ics_list, naming_list

def add_comparables(config,workflow):
    explain("what are comparables?...")
    config[workflow]['COMPARABLE'] = {}
    while True:
        a = conversation("do you want to add a comparable? enter 'y' or 'n'",None)
        if a=='y':
            name=conversation("enter the name",None)
            print("to remember, these are your ICS-keys:\n")
            for id in config[workflow].keys():
                if id == "TOOLS" or id == "COMPARABLE":
                    continue
                print(f"\t> '{id}'")
                for condition in config[workflow][id].keys():
                    print(f"\t> '{id}:{condition}'")
                    for setting in config[workflow][id][condition].keys():
                        print(f"\t> '{id}:{condition}:{setting}'")
            print("\n")
            group1=conversation("enter all ICS-keys for GROUP 1. ('id:condition:setting' or 'id:condition' or 'id')",None)
            group2=conversation("enter all ICS-keys for GROUP 2",None)
            config[workflow]['COMPARABLE'].update({name:[group1.split(','),group2.split(',')]})
        if a=='n':
            break
        else:
            print("enter 'y' or 'n'!")


def add_tools(config,workflow):
    tools_dict=config[workflow]['TOOLS']
    tools = conversation(f"choose from these tools: {list(tools_dict.keys())}",None)
    config[workflow]['TOOLS'] = {}
    for tool in tools.split(','):
        config[workflow]['TOOLS'].update({tool:tools_dict[tool]})


####################
####    MAIN    ####
####################

if __name__ == '__main__':

    makelogdir('LOGS')
    logid = scriptname+'.main: '
    # log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level=knownargs.loglevel)
    # log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog

    MIN_PYTHON = (3,7)
    if sys.version_info < MIN_PYTHON:
        log.error("This script requires Python version >= 3.7")
        sys.exit("This script requires Python version >= 3.7")
    # log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')

    print("\n\n","*"*30," SNAKES GUIDE ","*"* 30,"\n")

    try:
        args=parseargs()
        if args.manualmode:
            manual=True
            print("starting in manualmode\n\n")
        else:
            manual=False
            print("starting in quickmode\n\n")

        explain("Welcome...")
        start = conversation("enter 'append' for expanding an existing configfile, or enter 'new' for a new project",None)

        while True:

            if 'append' in start:
                configfile = conversation("enter the name of the existing file",None)
                config = load_configfile(os.path.abspath(configfile))
                print("\n")
                for key in get_path(config['SAMPLES']):
                    print("\t> ",key)
                print("\n")
                appending = conversation("this are all ICS i found in the file, would you like to add more? \nwrite 'y' for yes or 'n' for no", None)

                ics=""
                if 'y' in appending:
                    IDs = conversation("enter all ID's you want to work on",None)
                    if len(IDs.split(","))>1:
                        explain("okay, you specified mulitiple ID's, we will go through all the settings several times")
                    ICS=[]
                    for id in IDs.split(","):
                        conditions = conversation(f"now tell me all conditions for {id}",None)
                        for condi in conditions.split(","):
                            settings = conversation(f"and now all settings for {id}:{condi}",None)
                            for setting in settings.split(","):
                                ICS.append(f"{id}:{condi}:{setting}")
                    ics=",".join(ICS)

                print("okay, let's have a look to the workflows and postprocessing analyses")
                workflows = conversation("this are the current workflows, would you like to add more? Possible are MAPPING, TRIMMING, QC",config['WORKFLOWS'])
                postprocess = conversation("this are the current postprocess analyses, would you like to add more? Possible are DE, DEU, DAS",config['POSTPROCESSING'])

                annotation=""

                newconf= create_json_config(
                configfile=configfile,
                append="APPEND",
                template="snakes/configs/template_2.json",
                preprocess="",
                workflows=workflows,
                postprocess=postprocess,
                ics=ics,
                refdir=config['REFERENCE'],
                binaries=config['BINS'],
                procs=config['MAXTHREADS'],
                genomemap="",
                genomes="",
                genomeext="",
                sequencing="",
                annotation=annotation
                )
                break

            if 'new' in start:
                name = conversation("Please type the name of your Project, it will also be the name of the CONFIGFILE", None)
                configfile = f"config_{name}.json"

                explain("For each id to work on you can define one or multiple conditions and settings that will be used for the analysis). The ICS also sets the file structure to follow for the FASTQ directory, where the ID is the first level and the Condition the second. Setting is used by RunSnakemake.py to enable processing of the same samples under different settings like mapping tools, trimming tools and later also postprocessing tools or commandline options for these tools.")
                IDs = conversation("enter all ID's you want to work on",None)
                ICS=[]
                for id in IDs.split(","):
                    conditions = conversation(f"now tell me all conditions for {id}",None)
                    for condi in conditions.split(","):
                        settings = conversation(f"and now all settings for {id}:{condi}",None)
                        for setting in settings.split(","):
                            ICS.append(f"{id}:{condi}:{setting}")

                explain("Genomes, references...")
                fasta_dict={}
                organisms = conversation("enter all organisms you have in your analysis",None)
                for organism in organisms.split(","):
                    fasta_dict.update({organism:conversation(f"enter the FASTA.gz file appending to {organism}", None)})

                explain("workflow explanation...")
                workflows = conversation("which WORKFLOWS would you like to run? Possible are MAPPING, TRIMMING, QC", None)
                explain("postprocess explanation...")
                postprocess = conversation("which POSTPROCESS ANALYSIS would you like to run? Possible are DE, DEU, DAS", None)
                procs= conversation("enter the Maximum number of parallel processes to start snakemake with", None)

                annotation=""

                newconf = create_json_config(
                configfile = configfile,
                append="",
                template="snakes/configs/template_2.json",
                preprocess="",
                workflows=workflows,
                postprocess=postprocess,
                ics=",".join(ICS),
                refdir="GENOMES",
                binaries="snakes/scripts",
                procs=procs,
                genomemap=",".join([id+":organism" for id in IDs.split(",")]),
                genomes=",".join([k+":"+v for k,v in fasta_dict.items()]),
                genomeext="",
                sequencing="",
                annotation=annotation
                )
                break

            else:
                start = conversation("enter 'append' or 'new'",None)

        explain("we are now done with the basics. in the next steps we will set the details. for this we have to make the design of the experiment understandable for snakes")

        explain("genome extension...")
        rename(newconf,'NAME',
        ["enter the additional FASTA-extension"]
        )

        if len(newconf["GENOME"]) > 1:
            explain("source..")
            rename(newconf,'SOURCE',
            [f"enter the corresponding organism {list(newconf['GENOME'].keys())}"]
            )
        else:
            for id in newconf['SOURCE'].keys():
                for condition in newconf['SOURCE'][id].keys():
                    for setting in newconf['SOURCE'][id][condition].keys():
                        newconf['SOURCE'][id][condition][setting]=list(newconf['GENOME'].keys())[0]

        explain("sequencing types..")
        rename(newconf,'SEQUENCING',
        ["enter 'paired' or 'unpaired'"]
        )

        explain("samples...")
        for id in newconf['SAMPLES'].keys():
            for condition in newconf['SAMPLES'][id].keys():
                for setting in newconf['SAMPLES'][id][condition].keys():
                    samples = conversation(f"enter the sample-names of {id}:{condition}:{setting}",None)
                    newconf['SAMPLES'][id][condition][setting]=samples.split(",")

        explain("workflows...")

        if 'QC' in newconf['WORKFLOWS']:
            explain("QC...")
            rename(newconf,'QC',
            ["QC-options, would you like to change them?",'OPTIONS',0]
            )

        if 'MAPPING' in newconf['WORKFLOWS']:
            explain("MAPPING...")
            rename(newconf,'MAPPING',
            ["enter ANNOTATION file",'ANNOTATION'],
            ["INDEXING, would you like to change them?",'OPTIONS',0],
            ["MAPPER, would you like to change them?",'OPTIONS',1]
            )

        if 'TRIMMING' in newconf['WORKFLOWS']:
            explain("TRIMING...")
            rename(newconf,'TRIMMING',
            ["TRIMMING-options, would you like to change them?",'OPTIONS',0]
            )

        explain("now we come to the analysis settings..")

        if 'DE' in newconf['POSTPROCESSING']:
            explain("DE...")
            add_tools(newconf,'DE')
            add_comparables(newconf,'DE')
            rename(newconf,'DE',
            ["enter ANNOTATION file",'ANNOTATION'],
            ["Counting-settings for DE-Analysis, would you like to change them?",'OPTIONS',0]
            )

        if 'DEU' in newconf['POSTPROCESSING']:
            explain("DEU...")
            add_tools(newconf,'DEU')
            add_comparables(newconf,'DEU')
            rename(newconf,'DEU',
            ["enter ANNOTATION file",'ANNOTATION'],
            ["Counting-settings for DEU-Analysis, would you like to change them?",'OPTIONS',0]
            )

        if 'DAS' in newconf['POSTPROCESSING']:
            explain("DAS...")
            add_tools(newconf,'DAS')
            add_comparables(newconf,'DAS')
            rename(newconf,'DAS',
            ["enter ANNOTATION file",'ANNOTATION'],
            ["Counting-settings for DAS-Analysis, would you like to change them?",'OPTIONS',0],
            ["Analysis-options for diego for DAS-Analysis, would you like to change them?",'OPTIONS', 1]
            )

        print_json(newconf,configfile,annotation)

        print(f"ALL RIGHT! you will find now {configfile} for starting your analysis ")

        explain("How to start snakes with the configfile: ...")

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# Configurator.py ends here
