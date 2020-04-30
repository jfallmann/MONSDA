#!/usr/bin/env python3
# RunSnakemake.py ---
#
# Filename: RunSnakemake.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon Feb 10 08:09:48 2020 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Apr 30 22:15:13 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 947
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
#import snakemake
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re
#import logging
min_version("5.8.2")

from lib.Collection import *
from lib.Logger import *
scriptname=os.path.basename(__file__).replace('.py','')

def parseargs():
    parser = argparse.ArgumentParser(description='Wrapper around snakemake to run config based jobs automatically')
    parser.add_argument("-c", "--configfile", type=str, help='Configuration json to read')
    parser.add_argument("-g", "--debug-dag", action="store_true", help='Should the debug-dag be printed')
    parser.add_argument("-f", "--filegraph", action="store_true", help='Should the filegraph be printed')
    parser.add_argument("-d", "--directory", type=str, default='', help='Directory to work in')
    parser.add_argument("-u", "--use-conda", action="store_true", default=True, help='Should conda be used')
    parser.add_argument("-l", "--unlock", action="store_true", help='If directory is locked you can unlock before processing')
    parser.add_argument("-j", "--procs", type=int, default=1, help='Number of parallel processed to start snakemake with, capped by MAXTHREADS in config!')
    parser.add_argument("-s", "--skeleton", action="store_true", help='Just create the minimal directory hierarchy as needed')
    parser.add_argument("-v", "--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

def run_snakemake (configfile, debugdag, filegraph, workdir, useconda, procs, skeleton, loglevel, unlock=None, optionalargs=None):
    try:
        logid = scriptname+'.run_snakemake: '
        if skeleton:
            for subdir in ['SubSnakes', 'RAW', 'GENOMES', 'FASTQ', 'LOGS', 'TMP']:  # Add RAW for nanopore preprocessing
                makeoutdir(subdir)
            sys.exit('Skeleton directories created, please add files and rerun without --skeleton option')
        else:
            for subdir in ['SubSnakes', 'LOGS', 'TMP']:  # Add RAW for nanopore preprocessing
                makeoutdir(subdir)

        subdir = 'SubSnakes'
        config = load_configfile(configfile)
        argslist = list()
        if useconda:
            argslist.append("--use-conda")
        else:
            log.warning(logid+'You are not making use of conda, be aware that this will most likely not work for the workflows provided in this repository! To change append the --use-conda option to your commandline call. Tou can also preinstall all conda environments appending the --use-conda and the --create-envs-only arguments.')
        if debugdag:
            argslist.append("--debug-dag")
        if filegraph:
            argslist.append("--filegraph|dot|display")
        if optionalargs and len(optionalargs) > 0:
            log.debug(logid+'OPTIONALARGS: '+str(optionalargs))
            argslist.extend(optionalargs)
            if '--profile' in optionalargs and 'snakes/slurm' in optionalargs:
                makeoutdir('LOGS/SLURM')

        threads = min(int(config['MAXTHREADS']), procs) if 'MAXTHREADS' in config else procs

        if unlock:
            log.info(logid+'Unlocking directory')
            jobtorun = 'snakemake --unlock -j {t} -s {s} --configfile {c}'.format(t=threads, s=os.path.abspath(os.path.join('snakes','workflows','header.smk')), c=configfile)
            log.info(logid+'RUNNING '+str(jobtorun))
            job = runjob(jobtorun)
            log.debug(logid+'JOB CODE '+str(job))

        preprocess = subworkflows = postprocess = None

        if 'PREPROCESSING' in config:
            preprocess = config['PREPROCESSING'].split(',') # we keep this separate because not all preprocessing steps need extra configuration
            if len(preprocess) == 0 or preprocess[0] == '':
                preprocess = None
        if 'WORKFLOWS' in config:
            subworkflows = config['WORKFLOWS'].split(',')
            if len(subworkflows) == 0 or subworkflows[0] == '':
                subworkflows = None
        if 'POSTPROCESSING' in config:
            postprocess = config['POSTPROCESSING'].split(',') # we keep this separate because not all postprocessing steps need extra configuration
            if len(postprocess) == 0 or postprocess[0] == '':
                postprocess = None

        if preprocess:
            try:
                all([config[x] or x == '' for x in preprocess])
            except KeyError:
                log.warning(logid+'Not all required preprocessing steps have configuration in the config file')

        if subworkflows:
            try:
                all([config[x] or x == 'TRIMMING' or x == '' for x in subworkflows])
            except KeyError:
                log.warning(logid+'Not all required subworkflows have configuration in the config file')

        if postprocess:
            try:
                all([config[x] or x == '' for x in postprocess])
            except KeyError:
                log.warning(logid+'Not all required postprocessing steps have configuration in the config file')

        log.debug(logid+'WORKFLOWS: '+str([preprocess,subworkflows,postprocess]))

        '''
        Fix conda path if needed
        '''
        condapath=re.compile(r'conda:\s+"')
        logfix=re.compile(r'loglevel="INFO"')
        '''
        START TO PROCESS
        IF WE NEED TO DOWNLOAD FILES WE DO THIS NOW
        '''

        if preprocess and 'RAW' in preprocess:
            if 'RAW' not in config:
                log.error(logid+'No configuration with key \'RAW\' for file download found. Nothing to do!')
            makeoutdir('FASTQ')
            makeoutdir('TMP')
            preprocess.remove('RAW')
            SAMPLES = download_samples(config)
            log.info(logid+'PRESAMPLES: '+str(SAMPLES))
            conditions = get_conditions(SAMPLES,config) #[x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
            log.info(logid+'PRECONDITIONS: '+str(conditions))
            for condition in conditions:
                subconf = NestedDefaultDict()
                subwork = 'RAW'
                listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                if listoftools is None:
                    log.warning(logid+'No entry fits condition '+str(condition)+' for preprocessing step '+str(subwork))
                    continue
                toolenv, toolbin = map(str,listoftools[0])
                subconf.update(listofconfigs[0])
                subname = toolenv+'.smk'
                smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subsnake.smk'])))
                if os.path.exists(smko):
                    os.rename(smko,smko+'.bak')
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                            line = re.sub(condapath, 'conda:  "../', line)
                            smkout.write(line)
                    smkout.write('\n\n')

                smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                    smkout.write('\n\n')

                confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json'])))
                if os.path.exists(confo):
                    os.rename(confo,confo+'.bak')
                with open(confo, 'a') as confout:
                    json.dump(subconf, confout)

                jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                log.info(logid+'RUNNING '+str(jobtorun))
                job = runjob(jobtorun)
                log.debug(logid+'JOB CODE '+str(job))

        '''
        ONCE FILES ARE DOWNLOAD WE CAN START PROCESSING
        '''

        SAMPLES = get_samples(config)
        log.info(logid+'SAMPLES: '+str(SAMPLES))
        conditions = get_conditions(SAMPLES,config) #[x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
        log.info(logid+'CONDITIONS: '+str(conditions))

        if preprocess:
            log.info(logid+'STARTING PREPROCESSING')
            if 'QC' in preprocess and 'QC' in config:
                makeoutdir('QC')
            for condition in conditions:
                for subwork in preprocess:
                    subconf = NestedDefaultDict()
                    log.debug(logid+'PREPROCESS: '+str(subwork)+' CONDITION: '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    log.debug(logid+str([listoftools,listofconfigs]))
                    if listoftools is None:
                        log.warning(logid+'No entry fits condition '+str(condition)+' for preprocessing step '+str(subwork))
                        continue

                    for i in range(0,len(listoftools)):

                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.smk'
                        log.debug(logid+'PREPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                                #smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                        if subwork == 'QC':
                            subname = toolenv+'_raw.smk'

                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                        confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')
                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                        log.info(logid+'RUNNING '+str(jobtorun))
                        job = runjob(jobtorun)
                        log.debug(logid+'JOB CODE '+str(job))

        else:
            log.warning(logid+'No preprocessing workflows defined! Continuing with workflows!')

        if subworkflows:
            log.info(logid+'STARTING PROCESSING')
            for condition in conditions:
                smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk'])))
                if os.path.exists(smko):
                    os.rename(smko,smko+'.bak')
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                            line = re.sub(condapath,'conda:  "../',line)
                            smkout.write(line)
                            #smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                    smkout.write('\n\n')

                if 'QC' in subworkflows and 'QC' in config:
                    makeoutdir('QC')
                    smkf = os.path.abspath(os.path.join('snakes','workflows','multiqc.smk'))
                    if 'MAPPING' in subworkflows:
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write('rule themall:\n\tinput:\texpand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config)),\n\t\texpand("QC/Multi/{condition}/multiqc_report.html",condition=os.path.join(samplecond(SAMPLES,config)[0]))\n\n')
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')
                    else:
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write('rule themall:\n\tinput:\texpand("QC/Multi/{condition}/multiqc_report.html",condition=os.path.join(samplecond(SAMPLES,config)[0]))\n\n')
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                if 'MAPPING' in subworkflows and 'TRIMMING' not in subworkflows:
                    log.info(logid+'Simulating read trimming as trimming is not part of the workflow!')
                    makeoutdir('TRIMMED_FASTQ')
                    smkf = os.path.abspath(os.path.join('snakes','workflows','simulatetrim.smk'))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')

                subconf = NestedDefaultDict()
                for subwork in subworkflows:
                    log.debug(logid+'PREPARING '+str(subwork)+' '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.smk'
                        log.debug(logid+'SUBWORKFLOW: '+str([subwork,toolenv,subname,condition, subsamples, subconf]))

                        if subwork == 'QC' and 'TRIMMING' in subworkflows and not 'MAPPING' in subworkflows:
                            subname = toolenv+'_trim.smk'

                        if subwork == 'QC' and not 'TRIMMING' in subworkflows and not 'MAPPING' in subworkflows:
                            subname = toolenv+'_raw.smk'

                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                if 'MAPPING' in subworkflows:
                    smkf = os.path.abspath(os.path.join('snakes','workflows','mapping.smk'))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')

                confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json'])))
                if os.path.exists(confo):
                    os.rename(confo,confo+'.bak')
                with open(confo, 'a') as confout:
                    json.dump(subconf, confout)

            for condition in conditions:
                log.info(logid+'Starting workflows for condition '+str(condition))
                jobtorun = 'snakemake -j {t} -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                log.info(logid+'RUNNING WORKFLOW '+str(jobtorun))
                job = runjob(jobtorun)
                log.debug(logid+'JOB CODE '+str(job))

        else:
            log.warning(logid+'No subworkflows defined! Nothing to do!')

        if postprocess:
            log.info(logid+'STARTING POSTPROCESSING WITH SAMPLES '+str(SAMPLES))

            if 'PEAKS' in config and 'PEAKS' in postprocess:
                CLIP = checkclip(SAMPLES, config)
                log.info(logid+'Running Peak finding for '+CLIP+' protocol')

            for condition in conditions:
                subconf = NestedDefaultDict()
                for subwork in postprocess:
                    if any(subwork == x for x in ['DE', 'DEU', 'DAS']):
                        continue
                    log.debug(logid+'POSTPROCESS: '+str(subwork)+' CONDITION: '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    log.debug(logid+str([listoftools,listofconfigs]))
                    if listoftools is None:
                        log.warning(logid+'No entry fits condition '+str(condition)+' for postprocessing step '+str(subwork))
                        continue

                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subname = toolenv+'.smk'
                        subsamples = list(set(sampleslong(subconf)))
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))
                        smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                                #smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                        confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')

                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                        log.info(logid+'RUNNING '+str(jobtorun))
                        job = runjob(jobtorun)
                        log.debug(logid+'JOB CODE '+str(job))

            #THIS SECTION IS FOR DE, DEU, DAS ANALYSIS, WE USE THE CONDITIONS TO MAKE PAIRWISE COMPARISONS
            for analysis in ['DE', 'DEU', 'DAS']:
                if analysis in config and analysis in postprocess:
                    log.info(logid+'STARTING '+analysis+' Analysis...')
                    subwork = analysis
                    subconf = NestedDefaultDict()
                    log.debug(logid+'SUBWORK: '+str(subwork)+' CONDITION: '+str(conditions))
                    #listoftoolscount, listofconfigscount = create_subworkflow(config, 'COUNTING', conditions) #Counting is now done on per analysis rule to increase freedom for user
                    listoftools, listofconfigs = create_subworkflow(config, subwork, conditions)

                    if listoftools is None:# or listoftoolscount is None:
                        log.error(logid+'No entry fits condition '+str(conditions)+' for postprocessing step '+str(subwork))

                    for key in config[subwork]['TOOLS']:
                        log.info(logid+'... with Tool: '+key)
                        toolenv = key
                        toolbin = config[subwork]['TOOLS'][key]
                        #countenv, countbin = map(str,listoftoolscount[0]) #Counting per analysis rule now
                        subconf = NestedDefaultDict()
                        for i in listofconfigs:
                            i[subwork+'ENV'] = toolenv
                            i[subwork+'BIN'] = toolbin
                            #i['COUNTBIN'] = 'featureCounts'#This is hard coded where needed for now
                            #i['COUNTENV'] = 'countreads'#This is hard coded where needed for now
                        for i in range(len(listoftools)):
                            subconf = merge_dicts(subconf,listofconfigs[i])

                        #for x in range(0,len(listofconfigscount)): ### muss hier auch noch gefiltert werden?
                        #    subconf = merge_dicts(subconf,listofconfigscount[x])
                        subname = toolenv+'.smk' if toolenv != 'edger' else toolenv+'_'+subwork+'.smk'
                        subsamples = sampleslong(subconf)
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix,'loglevel="'+loglevel+'"',line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                                #smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subsnake.smk']))), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                        confo = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')
                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join([subwork,toolenv,'subsnake.smk'])]))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join([subwork,toolenv,'subconfig.json'])]))),d=workdir,rest=' '.join(argslist))
                        log.info(logid+'RUNNING '+str(jobtorun))
                        job = runjob(jobtorun)
                        log.debug(logid+'JOB CODE '+str(job))

        else:
            log.warning(logid+'No postprocessing steps defined! Nothing to do!')

        log.info('Workflows executed without error!')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))


def runjob(jobtorun):
    try:
        logid = scriptname+'.runjob: '
        #return subprocess.run(jobtorun, shell=True, universal_newlines=True, capture_output=True)  # python >= 3.7
        job = subprocess.Popen(jobtorun, shell=True, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, close_fds=True)

        while True:
            output = str.join('',job.stdout.readlines()).rstrip()
            err = str.join('',job.stderr.readlines()).rstrip()
            if output == '' and err == '' and job.poll() is not None:
                break
            if output and output != '':
                log.info(logid+str(output))
                if any(x in output for x in ['ERROR','Error','error','Exception']) and not 'Workflow finished' in output:
                    log.error(logid+'STOPPING: '+str(output))
                    job.kill()
                    sys.exit()
            if err and err != '':
                if not 'Workflow finished' in err and not 'Nothing to be done' in err and any(x in err for x in ['ERROR','Error','error','Exception']):
                    log.error(logid+'STOPPING: '+str(err))
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid+str(err))
            if job.poll() is not None:
                break

        if job.returncode == 0:
            output = str.join('',job.stdout.readlines()).rstrip()
            if output and output != '':
                log.info(logid+'JOB FINISHED: '+output)
            return job.poll()
        else:
            output = str.join('',job.stdout.readlines()).rstrip()
            err = str.join('',job.stderr.readlines()).rstrip()
            if err and err != '' or output and output != '':
                log.error(logid+'ERROR: '+errput+output)
            job.kill()
            sys.exit('ERROR SIGNAL: '+str(job.returncode))

    except Exception as err:
        job.kill()
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))
        sys.exit(''.join(tbe.format()))

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
        if not os.path.isfile(os.path.abspath('LOGS/'+scriptname+'.log')):
            open('LOGS/'+scriptname+'.log','a').close()
        log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level=knownargs.loglevel)
        log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=knownargs.loglevel)

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')
        log.debug(logid+str(log.handlers))

        run_snakemake(knownargs.configfile, knownargs.debug_dag, knownargs.filegraph, knownargs.directory, knownargs.use_conda, knownargs.procs, knownargs.skeleton, knownargs.loglevel, knownargs.unlock, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# RunSnakemake.py ends here
