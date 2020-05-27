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
import random
import logging

min_version("5.8.2")

from lib.Collection import *
from lib.Logger import *
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Helper to create initial config file used for workflow processing')
    parser.add_argument("-q", "--quickmode", action="store_true", help='If set configuration without explanation')

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
    oldics = []
    oldtodos = ['NAME','SOURCE','SEQUENCING','SAMPLES']

    todos = ','.join([x for x in [preprocess,workflows,postprocess] if x != '' ]).split(',')
    for x in todos:
        if x not in config and x != "":
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
            for id in oldconf['SAMPLES'].keys():
                for condition in oldconf['SAMPLES'][id].keys():
                    for setting in oldconf['SAMPLES'][id][condition].keys():
                        icslist.append(f"{id}:{condition}:{setting}")
                        oldics.append([[id],[condition],[setting]])
            for x in oldconf['PREPROCESSING'].split(","):
                oldtodos.append(x)
                todos.append(x)
            for x in oldconf['WORKFLOWS'].split(","):
                oldtodos.append(x)
                todos.append(x)
            for x in oldconf['POSTPROCESSING'].split(","):
                oldtodos.append(x)
                todos.append(x)
        if ics:
            for x in ics.split(','):
                if x not in icslist:
                    icslist.append(x)
    else:
        log.error(logid+'IdentifierConditionSetting (ics) not defined!')
        sys.exit()
    todos = [x for x in todos if x]

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
                    samplelist = get_samples_from_dir_2(id, condition, setting, newconf)
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
        else:
            newconf['PREPROCESSING'] = oldconf['PREPROCESSING']
        if workflows and workflows not in newconf['WORKFLOWS']:
            newconf['WORKFLOWS'] = str.join(',',list(set(str.join(',',[oldconf['WORKFLOWS'],workflows]).split(','))))
        else:
            newconf['WORKFLOWS'] = oldconf['WORKFLOWS']
        if postprocess and postprocess not in newconf['POSTPROCESSING']:
            newconf['POSTPROCESSING'] = str.join(',',list(set(str.join(',',[oldconf['POSTPROCESSING'],postprocess]).split(','))))
        else:
            newconf['POSTPROCESSING'] = oldconf['POSTPROCESSING']
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
        elif isinstance(oldconf['GENOME'], dict):
            newconf['GENOME'] = oldconf['GENOME']
        else:
            newconf['GENOME'] = str(oldconf['GENOME'])

        log.debug(logid+'GENOMEMAPCONF: '+str(newconf['GENOME']))

        for key in ['NAME','SOURCE','SEQUENCING','SAMPLES']:
            for id,condition,setting in [x.split(':') for x in icslist]:
                try:
                    newconf[key][id][condition][setting] = oldconf[key][id][condition][setting]
                except:
                    if key == 'SAMPLES':
                        samplelist = get_samples_from_dir_2(id, condition, setting, newconf)
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
                newconf[key][id][condition][setting].update(config[key]['id']['condition']['setting'])
                log.debug(logid+'TODO: '+str(key)+'\t'+str(config[key])+'\t'+str(newconf[key]))

            if 'id' in newconf[key]:
                newconf[key][id] = newconf[key].pop('id')
                newconf[key][id][condition] = newconf[key][id].pop('condition')
                newconf[key][id][condition][setting] = newconf[key][id][condition].pop('setting')

        if key=='DE' or key=="DEU" or key=='DAS':
            set_relations(newconf,key)

    input=json.dumps(newconf)
    flatconf=json.loads(input)

    return oldtodos, oldics, flatconf


def set_relations(config,key):
    for id in config['SAMPLES'].keys():
        for condition in config['SAMPLES'][id].keys():
            for setting in config['SAMPLES'][id][condition].keys():
                if config["SAMPLES"][id][condition][setting]:
                    relations=[]
                    ics = f"{id}:{condition}:{setting}"
                    type = config["SEQUENCING"][id][condition][setting]
                    for sample in config["SAMPLES"][id][condition][setting]:
                        relations.append((sample,ics,type))
                    config[key][id][condition][setting]['RELATIONS']=relations

# @check_run
def print_json(paramdict,ofn):
    with open(os.path.join(project,ofn),'w') as jsonout:
        print(json.dumps(paramdict,indent=4),file=jsonout)

def proof_input(proof=None):
    allowed_characters=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','1','2','3','4','5','6','7','8','9','0','(',')','_','-',',','.',':','/']
    while True:
        a = input(">>> ").strip().replace(" ","")
        # print("\n")
        if any(x not in allowed_characters for x in a):
            print("You used unallowed letters, try again")
            continue
        if proof is not None and proof != "only-numbers" and any(x not in proof for x in a.split(",")):
            print(f"available are only: {proof}")
            continue
        if proof=="only-numbers":
            try:
                float(a)
                return a
            except:
                print("please enter integer or float")
                continue
        else:
            return a

def conversation(question, origin, proof=None):
    print(question)
    default = copy.deepcopy(origin)
    if default is None:
        a = proof_input(proof)
        print("\n")
        return a
    else:
        print("\n")

        if isinstance(default, str):
            print("\t> ",default)
            print("\n")
            if proof:
                print("enter what should be added ")
            else:
                print("enter what should be added or enter 'n' to continue ")
            a = proof_input(proof)
            print("\n")
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
                    return default
                if a=='y':
                    print("okay, tell me your setting:")
                    if default:
                        for element in default:
                            print(element+":")
                            a = proof_input()
                            default[element]=a
                        return default
                    else:
                        while True:
                            print("enter the name of the argument or enter 'n' to quit")
                            name = proof_input()
                            if name=='n':
                                break
                            print("enter the value")
                            value = proof_input()
                            print("\n")
                            default.update({name:value})
                    return default

def list_generator(config,workflow):
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

def rename(config,job,oldics,oldtodos,proof,*args):
    ics_list,naming_list = list_generator(config,job)
    if job in oldtodos:
        ics_list = [x for x in ics_list if x not in oldics]
        naming_list = [x for x in naming_list if x not in oldics]
    if not ics_list:
        return
    for naming in naming_list:
        switch=False
        for ics in ics_list:
                if naming==ics:
                    switch=True
                    print(f"okay, look at the settings of {ics[0][0]}:{ics[1][0]}:{ics[2][0]}")
                    for call in args:
                        if len(call)==1:
                            text  =call[0]
                            config[job][ics[0][0]][ics[1][0]][ics[2][0]] = conversation(text,None,proof)
                        if len(call)==2:
                            text  =call[0]
                            key   =call[1]
                            if 'ANNOTATION' in key:
                                org = config['SOURCE'][ics[0][0]][ics[1][0]][ics[2][0]]
                                files = os.listdir(os.path.join(project,'GENOMES',org))
                                if len(files) > 0:
                                    print(f"Following annotation-files are available for {org}:")
                                    print("\n")
                                    for file in files:
                                        print("\t> ",os.path.basename(file))
                                    print("\n")
                                    while True:
                                        anno = conversation(text,None,proof)
                                        if anno not in files:
                                            a = conversation(f"WARNING: {anno} is not yet in /GENOMES/{org} \nEnter 'y' for enter it anyway or enter 'n' to try again.",None)
                                            if 'n' in a:
                                                continue
                                            if 'y' in a:
                                                config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = anno
                                                break
                                        config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = anno
                                        break
                                else:
                                    print(f"There are no annotation files available for {org}. you can still specify \na file, but then remember to copy it to the directory")
                                    config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = conversation(text,None,proof)
                            else:
                                config[job][ics[0][0]][ics[1][0]][ics[2][0]][key] = conversation(text,None,proof)
                        if len(call)==3:
                            text  =call[0]
                            key   =call[1]
                            number=call[2]
                            config[job][ics[0][0]][ics[1][0]][ics[2][0]][key][number] = conversation(text, config[job][naming[0][0]][naming[1][0]][naming[2][0]][key][number],proof)
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



def get_path(input_dict):
    for key, value in input_dict.items():
        if isinstance(value, dict):
            for subkey in get_path(value):
                yield key + ':' + subkey
        else:
            yield key

def add_comparables(config,workflow):
    config[workflow]['COMPARABLE'] = {}
    while True:
        a = conversation("do you want to add a COMPARABLE? enter 'y' or 'n'",None)
        if a=='y':
            name=conversation("enter the name",None)
            ics_list=[]
            print("to remember, these are your ICS-keys:\n")
            for id in config[workflow].keys():
                if id == "TOOLS" or id == "COMPARABLE":
                    continue
                ics = f"{id}"
                ics_list.append(ics)
                print(f"\t> {id}")
                for condition in config[workflow][id].keys():
                    ics = f"{id}:{condition}"
                    ics_list.append(ics)
                    print(f"\t> {id}:{condition}")
                    for setting in config[workflow][id][condition].keys():
                        ics = f"{id}:{condition}:{setting}"
                        ics_list.append(ics)
                        print(f"\t> {ics}")
            print("\n")
            group1=conversation("enter all ICS-keys for GROUP 1.",None,ics_list)
            group2=conversation("enter all ICS-keys for GROUP 2",None,ics_list)
            config[workflow]['COMPARABLE'].update({name:[group1.split(','),group2.split(',')]})
            print("\n")
        elif a=='n':
            break
        else:
            print("enter 'y' or 'n'!")

def add_tools(config,workflow):
    tools_dict=config[workflow]['TOOLS']
    tools = conversation(f"choose from these tools: {list(tools_dict.keys())}",None)
    config[workflow]['TOOLS'] = {}
    for tool in tools.split(','):
        config[workflow]['TOOLS'].update({tool:tools_dict[tool]})

def explain(filename):
    if not quick:
        path = os.path.join("docs/guide",filename)
        with open(path, 'r') as f:
            text = f.read()
        print("\n")
        print(text,"\n")

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

    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    cwd=os.getcwd()

    print("\n\n\n","*"*33," NextSnakes GUIDE ","*"* 33,)
    print("","*"*86,"\n")

    try:
        args=parseargs()
        if args.quickmode:
            quick=True
            print("  > starting in quick-mode <\n")
        else:
            quick=False
            print("  > starting in explanation-mode <\n")

        explain("intro.txt")
        start = conversation("Enter 'append' for expanding an existing configfile or 'new' for a new project",None)

        while True:

            if 'append' in start:

                while True:
                    configfile = conversation("enter the name of the existing file",None)
                    if os.path.isfile(configfile):
                        break
                    else:
                        print("Sry, didn't find such a file...")
                config = load_configfile(os.path.abspath(configfile))
                print("\n")
                for key in get_path(config['SAMPLES']):
                    print("\t> ",key)
                print("\n")
                appending = conversation("this are all ICS i found in the file, would you like to add more? \nwrite 'y' for yes or 'n' for no", None,['y','n'])

                ICS=[]
                if 'y' in appending:
                    IDs = conversation("enter all ID's you want to work on",None)
                    for id in IDs.split(","):
                        conditions = conversation(f"now tell me all conditions for {id}",None)
                        for condi in conditions.split(","):
                            settings = conversation(f"and now all settings for {id}:{condi}",None)
                            for setting in settings.split(","):
                                ICS.append(f"{id}:{condi}:{setting}")

                print("okay, let's have a look to the workflows and postprocessing analyses")
                workflows = conversation("this are the current workflows, would you like to add more? Possible are MAPPING, TRIMMING, QC",config['WORKFLOWS'],["MAPPING","TRIMMING","QC",""])
                postprocess = conversation("this are the current postprocess analyses, would you like to add more? Possible are DE, DEU, DAS",config['POSTPROCESSING'],["DE","DEU","DAS",""])

                annotation=""

                oldtodos, oldics, newconf= create_json_config(
                configfile=configfile,
                append="APPEND",
                template="configs/template_2.json",
                preprocess="",
                workflows=workflows,
                postprocess=postprocess,
                ics=",".join(ICS),
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

                name = conversation("Now please type the name of your Project, it will be the name of the CONFIGFILE \nand possibly of your Project-folder", None)
                configfile = f"config_{name}.json"

                folder_content = os.listdir('../')
                if 'FASTQ' in folder_content or 'GENOMES' in folder_content:
                    set_folder=conversation("It looks like you already set up your project-folder. We would therefor skip setting it up now. Enter 'n' if you want to do that anyway.", None)
                else:
                    explain('projectfolder.txt')
                    while True:
                        path_to_project = conversation("So, where should the Project-folder be created? Enter the absolute path \nlike '/home/path/to/NextSnakesProjects' ",None)
                        if os.path.isdir(path_to_project):
                            project = os.path.join(path_to_project,name)
                            if os.path.isdir(project):
                                print(f"In the directory you entered, a folder with the name {name} already exist. \nMaybe you should think about what you really want to do. If you want to work on \nan existing Project, use the 'append' function, otherwise use a different \nProject-name. Ciao!")
                                exit()
                            os.mkdir(project)
                            os.mkdir(os.path.join(project,"FASTQ"))
                            os.mkdir(os.path.join(project,"GENOMES"))
                            os.symlink(cwd, os.path.join(project,'snakes'))
                            break
                        else:
                            print("Sry, couldn't find this directory ")

                explain("ics.txt")
                IDs = conversation("enter all ID's you want to work on",None)
                ICS=[]
                for id in IDs.split(","):
                    os.mkdir(os.path.join(project,"FASTQ",id))
                    conditions = conversation(f"now tell me all conditions for {id}",None)
                    for condi in conditions.split(","):
                        os.mkdir(os.path.join(project,"FASTQ",id,condi))
                        settings = conversation(f"and now all settings for {id}:{condi}",None)
                        for setting in settings.split(","):
                            os.mkdir(os.path.join(project,"FASTQ",id,condi,setting))
                            ICS.append(f"{id}:{condi}:{setting}")

                explain("workflows.txt")
                workflows = conversation("Enter which WORKFLOWS you would like to run. Possible are MAPPING, TRIMMING, QC", None,["MAPPING", "TRIMMING", "QC",""])
                postprocess = conversation("Which POSTPROCESS ANALYSIS would you like to run? Possible are DE, DEU, DAS", None, ["DE","DEU","DAS",""])

                print("Now you may have to be patient. We go through all settings, and you have to \n"
                    "assign the corresponding samples by their absolute directory path. The Guide will \n"
                    "symlink it to your 'FASTQ'-Folder.\n\n"
                    "For that, maybe it's helpful to open another termimal, navigate to the directory \n"
                    "your samples are stored and list them with 'readlink -f *.fastq.gz'\n")

                for ics in ICS:
                    while True:
                        sample=conversation(f"Enter one Sample-path of {ics} or enter 'n' to go to the next ics",None)
                        if sample=='n':
                            break
                        if os.path.isfile(sample):
                            try:
                                os.symlink(sample, os.path.join(project,'FASTQ',ics.split(':')[0],ics.split(':')[1],ics.split(':')[2],os.path.basename(sample)))
                            except:
                                print("hmm, an error occured.. maybe this file is already linked. try again!")
                        else:
                            print("Sry, couldn't find the file")

                explain("genomes.txt")
                fasta_dict={}
                organisms = conversation("enter all organisms you have in your analysis",None)
                for organism in organisms.split(","):
                    os.mkdir(os.path.join(project,"GENOMES",organism))
                    fasta_dict.update({organism:conversation(f"enter the basename(!) of the fa.gz file appending to {organism}", None)})
                    print("Now you can add Genome reference files to your GENOMES folder\n")
                    while True:
                        file=conversation(f"Enter one file-path you want to symlink to {organism} or enter 'n' to go to the next ics",None)
                        if file=='n':
                            break
                        if os.path.isfile(file):
                            try:
                                os.symlink(file, os.path.join(project,'GENOMES',organism,os.path.basename(file)))
                            except:
                                print("hmm, an error occured.. maybe this file is already linked. try again!")
                        else:
                            print("Sry, couldn't find the file")


                annotation=""

                oldtodos, oldics, newconf = create_json_config(
                configfile = configfile,
                append="",
                template="configs/template_2.json",
                preprocess="",
                workflows=workflows,
                postprocess=postprocess,
                ics=",".join(ICS),
                refdir="GENOMES",
                binaries="scripts",
                procs="10",
                genomemap=",".join([id+":organism" for id in IDs.split(",")]),
                genomes=",".join([k+":"+v for k,v in fasta_dict.items()]),
                genomeext="",
                sequencing="",
                annotation=annotation
                )
                break

            else:
                start = conversation("enter 'append' or 'new'",None)

        # print(newconf)

        if ICS:
            explain("genomes2.txt")
            # SOURCE
            if len(newconf["GENOME"]) > 1:
                rename(newconf,'SOURCE',oldics,oldtodos,list(newconf['GENOME'].keys()),
                [f"enter the corresponding organism {list(newconf['GENOME'].keys())}"]
                )
            else:
                for id in newconf['SOURCE'].keys():
                    for condition in newconf['SOURCE'][id].keys():
                        for setting in newconf['SOURCE'][id][condition].keys():
                            newconf['SOURCE'][id][condition][setting]=list(newconf['GENOME'].keys())[0]

            # NAME
            rename(newconf,'NAME',oldics,oldtodos,None,
            ["enter the additional FASTA-extension"]
            )

            # SEQUENCING
            rename(newconf,'SEQUENCING',oldics,oldtodos,["paired", "unpaired"],
            ["enter 'paired' or 'unpaired'"]
            )

            # # SAMPLES
            # explain("samples...")
            # for id in newconf['SAMPLES'].keys():
            #     for condition in newconf['SAMPLES'][id].keys():
            #         for setting in newconf['SAMPLES'][id][condition].keys():
            #             if not newconf['SAMPLES'][id][condition][setting]:
            #                 samples = conversation(f"I couldn't found samples for {id}:{condition}:{setting}, please type the file names",None)
            #                 newconf['SAMPLES'][id][condition][setting]=samples.split(",")

        explain("workflows2.txt")

        if 'QC' in newconf['WORKFLOWS'].split(",") and ('QC' not in oldtodos or ICS):
            # explain("QC...")
            rename(newconf,'QC',oldics,oldtodos,None,
            ["QC-options, would you like to change them?",'OPTIONS',0]
            )

        if 'MAPPING' in newconf['WORKFLOWS'].split(",") and ('MAPPING' not in oldtodos or ICS):
            # explain("MAPPING...")
            rename(newconf,'MAPPING',oldics,oldtodos,None,
            ["enter ANNOTATION file",'ANNOTATION'],
            ["INDEXING, would you like to change them?",'OPTIONS',0],
            ["MAPPER, would you like to change them?",'OPTIONS',1]
            )

        if 'TRIMMING' in newconf['WORKFLOWS'].split(",") and ('TRIMMING' not in oldtodos or ICS):
            # explain("TRIMING...")
            rename(newconf,'TRIMMING',oldics,oldtodos,None,
            ["TRIMMING-options, would you like to change them?",'OPTIONS',0]
            )

        explain("analysis.txt")

        if 'DE' in newconf['POSTPROCESSING'].split(",") and ('DE' not in oldtodos or ICS):
            # explain("DE...")
            set_relations(newconf,'DE')
            add_tools(newconf,'DE')
            add_comparables(newconf,'DE')
            rename(newconf,'DE',oldics,oldtodos,None,
            ["enter ANNOTATION file (consider, that it has to be )",'ANNOTATION'],
            ["Counting-settings for DE-Analysis, would you like to change them?",'OPTIONS',0]
            )

        if 'DEU' in newconf['POSTPROCESSING'].split(",") and ('DEU' not in oldtodos or ICS):
            # explain("DEU...")
            set_relations(newconf,'DEU')
            add_tools(newconf,'DEU')
            add_comparables(newconf,'DEU')
            rename(newconf,'DEU',oldics,oldtodos,None,
            ["enter ANNOTATION file",'ANNOTATION'],
            ["Counting-settings for DEU-Analysis, would you like to change them?",'OPTIONS',0]
            )

        if 'DAS' in newconf['POSTPROCESSING'].split(",") and ('DAS' not in oldtodos or ICS):
            # explain("DAS...")
            set_relations(newconf,'DAS')
            add_tools(newconf,'DAS')
            add_comparables(newconf,'DAS')
            rename(newconf,'DAS',oldics,oldtodos,None,
            ["enter ANNOTATION file",'ANNOTATION'],
            ["Counting-settings for DAS-Analysis, would you like to change them?",'OPTIONS',0],
            ["Analysis-options for diego for DAS-Analysis, would you like to change them?",'OPTIONS', 1]
            )

        explain("processes.txt")
        procs= conversation("Enter the Maximum number of parallel processes to start snakemake with", None, "only-numbers")

        print_json(newconf,configfile)

        explain("runsnakemake.txt")
        print(f"\nYou created \n\n\t> {configfile}\n")
        print(f"start RunSnakemake with:\n\n\t> python3 snakes/RunSnakemake.py -c {configfile} --directory ${{PWD}}\n")
        print("\nCIAO!")


    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# Configurator.py ends here
