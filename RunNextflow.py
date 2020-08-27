#!/usr/bin/env python3
# RunNextflow.py ---
#
# Filename: RunNextflow.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Mon May 18 08:09:48 2020 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Aug 27 10:09:01 2020 (+0200)
#           By: Joerg Fallmann
#     Update #: 1411
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
import argparse
import subprocess
import re

scriptname=os.path.basename(__file__).replace('.py','')

from lib.Logger import *
#Logging
import datetime
makelogdir('LOGS')
if not os.path.isfile(os.path.abspath('LOGS/'+scriptname+'.log')):
    open('LOGS/'+scriptname+'.log','a').close()
else:
    ts = str(datetime.datetime.fromtimestamp(os.path.getmtime(os.path.abspath('LOGS/'+scriptname+'.log'))).strftime("%Y%m%d_%H_%M_%S"))
    shutil.copy2('LOGS/'+scriptname+'.log','LOGS/'+scriptname+'_'+ts+'.log')

log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M')
log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M')

#import Collection
from lib.Collection import *

def parseargs():
    parser = argparse.ArgumentParser(description='Wrapper around nextflow to run config based jobs automatically')
    parser.add_argument("-c", "--configfile", type=str, help='Configuration json to read')
    parser.add_argument("-d", "--directory", type=str, default='${PWD}', help='Directory to work in')
    parser.add_argument("-j", "--procs", type=int, default=1, help='Number of parallel processed to start nextflow with, capped by MAXTHREADS in config!')
    parser.add_argument("--clean", action="store_true", help='Cleanup workdir, append -n to see list of files to clean or -f to actually remove those files')
    parser.add_argument("-s", "--skeleton", action="store_true", help='Just create the minimal directory hierarchy as needed')
    parser.add_argument("-v", "--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

def run_nextflow (configfile, workdir, procs, skeleton, loglevel, clean=None, optionalargs=None):
    try:
        logid = scriptname+'.run_nextflow: '
        argslist = list()
        if optionalargs and len(optionalargs) > 0:
            log.debug(logid+'OPTIONALARGS: '+str(optionalargs))
            argslist.extend(optionalargs)

        if clean:
            log.info(logid+'Cleaning working directory')
            jobtorun = 'nextflow clean {a}'.format(a=' '.join(argslist))
            log.info(logid+'CLEANUP '+str(jobtorun))
            job = runjob(jobtorun)
            log.debug(logid+'JOB CODE '+str(job))
            sys.exit()

        config = load_configfile(configfile)
        refdir = config['REFERENCES']
        workdir = os.path.abspath(str.join(os.sep,[workdir,'NextFlowWork']))
        if skeleton:
            for subdir in ['SubSnakes', 'RAW', refdir, 'FASTQ', 'LOGS', 'TMP']:  # Add RAW for nanopore preprocessing
                makeoutdir(subdir)
            sys.exit('Skeleton directories created, please add files and rerun without --skeleton option')
        else:
            for subdir in ['SubFlows', 'LOGS', 'TMP', workdir]:  # Add RAW for nanopore preprocessing
                makeoutdir(subdir)

        subdir = 'SubFlows'

        argslist = list()
        if optionalargs and len(optionalargs) > 0:
            log.debug(logid+'OPTIONALARGS: '+str(optionalargs))
            argslist.extend(optionalargs)
            if '--profile' in optionalargs and 'nextsnakes/slurm' in optionalargs:  # NEEDS REFIT FOR NEXTFLOW
                makeoutdir('LOGS/SLURM')

        threads = min(int(config['MAXTHREADS']), procs) if 'MAXTHREADS' in config else procs

        preprocess = subworkflows = postprocess = None

        #Define workflow stages
        pre = ['QC','RAW']
        sub = ['QC','MAPPING','TRIMMING']
        post = ['COUNTING','UCSC','PEAKS','DE','DEU','DAS','ANNOTATE']

        wfs = config['WORKFLOWS'].split(',')

        if 'WORKFLOWS' in config:
            subworkflows = [x for x in wfs if x in sub]
            if len(subworkflows) == 0 or subworkflows[0] == '':
                subworkflows = None
            for x in subworkflows:
                wfs.remove(x)
            preprocess = [x for x in wfs if x in pre]
            if len(preprocess) == 0 or preprocess[0] == '':
                preprocess = None
            postprocess = [x for x in wfs if x in post]
            if len(postprocess) == 0 or postprocess[0] == '':
                postprocess = None
        else:
            log.error('NO WORKFLOWS DEFINED, NOTHING TO DO!')
            sys.exit()

        log.info('PRE: '+str(preprocess)+' SUB: '+str(subworkflows) + ' POST: '+str(postprocess))

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

        condapath=re.compile(r'conda\s+\"')

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
            conditions = get_conditions(SAMPLES,config)
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
                subsamples = list(set(sampleslong(subconf)))
                subname = toolenv+'.nf'
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.nf'))
                smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subflow.nf'])))
                if os.path.exists(smko):
                    os.rename(smko,smko+'.bak')
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            smkout.write(line)
                    smkout.write('\n\n')

                smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())
                    smkout.write('\n\n')

                #workflow merger
                with open(smko, 'a') as smkout:
                    smkout.write('\n\n'+'workflow {\n    main:\n')
                    for w in ['RAW']:
                        if w in flowlist:
                            smkout.write('\n    '+w+'()')
                        smkout.write('\n}')

                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.nf'))
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())

                confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json'])))
                if os.path.exists(confo):
                    os.rename(confo,confo+'.bak')
                with open(confo, 'a') as confout:
                    json.dump(subconf, confout)

                params = nf_fetch_params(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subconfig.json']))))
                toolparams = nf_tool_params(subsamples[0], None, subconf, subwork, toolenv, toolbin)

                jobtorun = 'nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subflow.nf']))), d=workdir, rest=' '.join(argslist), p=' '.join("--{!s} {!s}".format(key,val) for (key,val) in params.items()), j=toolparams, c = '--CONDITION '+str.join(os.sep,condition))

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
                log.debug(logid+'Working on condition: '+str(condition))
                for subwork in preprocess:
                    subconf = NestedDefaultDict()
                    log.debug(logid+'PREPROCESS: '+str(subwork)+' CONDITION: '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    log.debug(logid+str([listoftools,listofconfigs]))
                    if listoftools is None:
                        log.warning(logid+'No entry fits condition '+str(condition)+' for preprocessing step '+str(subwork))
                        continue

                    for i in range(0,len(listoftools)):
                        flowlist = list()
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.nf'
                        log.debug(logid+'PREPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.nf'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subflow.nf'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    #line = re.sub(condapath,'conda  \"../',line)
                                    smkout.write(line)
                            smkout.write('\n\n')

                        if subwork == 'QC':
                            subname = toolenv+'_raw.nf'
                            flowlist.append('QC_RAW')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                                smkout.write(smk.read())
                            smkout.write('\n\n')

                        if subwork == 'QC':
                            smkf = os.path.abspath(os.path.join('nextsnakes','workflows','multiqc.nf'))
                            with open(smko, 'a') as smkout:
                                with open(smkf,'r') as smk:
                                    smkout.write(smk.read())
                                    #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                                smkout.write('\n\n')
                            flowlist.append('MULTIQC')

                        #workflow merger
                        with open(smko, 'a') as smkout:
                            smkout.write('\n\n'+'workflow {\n    main:\n')
                            for w in ['QC_RAW']:
                                if w in flowlist:
                                    smkout.write('\n    '+w+'()')
                                smkout.write('\n}')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.nf'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())

                        confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')
                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        params = nf_fetch_params(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subconfig.json']))))
                        toolparams = nf_tool_params(subsamples[0], None, subconf, subwork, toolenv, toolbin)

                        jobtorun = 'nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subflow.nf']))), d=workdir, rest=' '.join(argslist), p=' '.join("--{!s} {!s}".format(key,val) for (key,val) in params.items()), j=toolparams, c = '--CONDITION '+str.join(os.sep,condition))

                        log.info(logid+'RUNNING '+str(jobtorun))
                        job = runjob(jobtorun)
                        log.debug(logid+'JOB CODE '+str(job))

        else:
            log.warning(logid+'No preprocessing workflows defined! Continuing with workflows!')

        '''
        END OF PREPROCESSING, START OF PROCESSING
        '''

        if subworkflows:
            log.info(logid+'STARTING PROCESSING')
            for condition in conditions:
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.nf'))
                smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subflow.nf'])))
                flowlist = list()
                if os.path.exists(smko):
                    os.rename(smko,smko+'.bak')
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            #line = re.sub(condapath,'conda  \"../',line)
                            smkout.write(line)
                    smkout.write('\n\n')

                if 'QC' in subworkflows and 'QC' in config:
                    makeoutdir('QC')

                if 'MAPPING' in subworkflows and 'QC' not in subworkflows:
                    log.info(logid+'Mapping without QC!')

                if 'MAPPING' in subworkflows and 'TRIMMING' not in subworkflows:
                    log.info(logid+'Simulating read trimming as trimming is not part of the workflow!')
                    makeoutdir('TRIMMED_FASTQ')
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','simulatetrim.nf'))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                            smkout.write(smk.read())
                        smkout.write('\n\n')

                if 'TRIMMING' in subworkflows and 'QC' not in subworkflows and 'MAPPING' not in subworkflows:
                    log.info(logid+'Trimming without QC!')

                subconf = NestedDefaultDict()
                subnames = list()
                for subwork in subworkflows:
                    log.debug(logid+'PREPARING '+str(subwork)+' '+str(condition))
                    if subwork != 'QC':
                        flowlist.append(subwork)
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.nf'
                        if subwork != 'QC':
                            subnames.append(subname)
                        log.debug(logid+'SUBWORKFLOW: '+str([subwork,toolenv,subname,condition, subsamples, subconf]))

                        if subwork == 'QC':
                            if 'TRIMMING' in subworkflows:
                                subnames.extend([toolenv+'_raw.nf', toolenv+'_trim.nf'])
                                flowlist.extend(['QC_RAW','QC_TRIMMING'])
                            else:
                                subnames.append(toolenv+'_raw.nf')
                                flowlist.append('QC_RAW')
                            if 'MAPPING' in subworkflows:
                                subnames.append(toolenv+'.nf')
                                flowlist.append('QC_MAPPING')

                if 'MAPPING' in subworkflows:
                    subnames.append('mapping.nf')
                if 'QC' in subworkflows:
                    subnames.append('multiqc.nf')
                    flowlist.append('MULTIQC')

                for subname in subnames:
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                            smkout.write(smk.read())
                        smkout.write('\n\n')

                #workflow merger
                log.debug('FLOWLIST: '+str(flowlist))
                paired = checkpaired([SAMPLES[0]],config)
                with open(smko, 'a') as smkout:
                    smkout.write('\n\n'+'workflow {\n')
                    for w in ['QC_RAW','TRIMMING','QC_TRIMMING','MAPPING','QC_MAPPING','MULTIQC']:
                        if w in flowlist:
                            if w ==  'QC_TRIMMING':
                                smkout.write(' '*4+w+'(TRIMMING.out.trimmed)\n')
                            elif w == 'MAPPING':
                                smkout.write(' '*4+w+'(TRIMMING.out.trimmed)\n')
                                smkout.write(' '*4+'POSTMAPPING(MAPPING.out.mapped)\n')
                            elif w == 'QC_MAPPING':
                                smkout.write(' '*4+w+'(POSTMAPPING.out.postmapuni)\n')
                            elif w ==  'MULTIQC':
                                if 'MAPPING' in flowlist:
                                    smkout.write(' '*4+w+'(QC_MAPPING.out.qc)\n')
                                elif 'TRIMMING' in flowlist:
                                    smkout.write(' '*4+w+'(QC_TRIMMING.out.qc)\n')
                                else:
                                    smkout.write(' '*4+w+'(QC_RAW.out.qc)\n')
                            else:
                                smkout.write(' '*4+w+'(dummy)\n')
                    smkout.write('}\n\n')

                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.nf'))
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())

                confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json'])))
                if os.path.exists(confo):
                    os.rename(confo,confo+'.bak')
                with open(confo, 'a') as confout:
                    json.dump(subconf, confout)

            for condition in conditions:
                log.info(logid+'Starting workflows for condition '+str(condition))

                params = nf_fetch_params(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))))
                toolparams = nf_tool_params(subsamples[0], None, config, subwork, toolenv, toolbin, subworkflows, condition)

                jobtorun = 'nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subflow.nf']))), d=workdir, rest=' '.join(argslist), p=' '.join("--{!s} {!s}".format(key,val) for (key,val) in params.items()), j=toolparams, c = '--CONDITION '+str.join(os.sep,condition))

                log.info(logid+'RUNNING '+str(jobtorun))
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
                        subname = toolenv+'.nf'
                        subsamples = list(set(sampleslong(subconf)))
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.nf'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subflow.nf'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    #line = re.sub(condapath,'conda  \"../',line)
                                    smkout.write(line)
                                #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                                smkout.write(smk.read())
                            smkout.write('\n\n')

                        #workflow merger
                        with open(smko, 'a') as smkout:
                            smkout.write('\n\n'+'workflow {\n    main:\n')
                            for w in ['QC_RAW','TRIMMING','QC_TRIMMING','MAPPING','QC_MAPPING','MULTIQC']:
                                if w in flowlist:
                                    smkout.write('\n    '.w)
                                smkout.write('\n}')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.nf'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())

                        confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')

                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        params = nf_fetch_params(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))))
                        toolparams = nf_tool_params(subsamples[0], None, subconf, subwork, condition)

                        jobtorun = 'nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subflow.nf']))), d=workdir, rest=' '.join(argslist), p=' '.join("--{!s} {!s}".format(key,val) for (key,val) in params.items()), j=toolparams, c = '--CONDITION '+str.join(os.sep,condition))

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
                        subname = toolenv+'.nf' if toolenv != 'edger' else toolenv+'_'+subwork+'.nf'
                        subsamples = sampleslong(subconf)
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.nf'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subflow.nf'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    #line = re.sub(condapath,'conda  \"../',line)
                                    smkout.write(line)
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subflow.nf']))), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                #smkout.write(re.sub(condapath,'conda  \"../',smk.read()))
                                smkout.write(smk.read())
                            smkout.write('\n')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.nf'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())

                        confo = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')
                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        params = nf_fetch_params(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))))
                        toolparams = nf_tool_params(subsamples[0], None, subconf, subwork, condition)

                        jobtorun = 'nextflow -log /dev/stderr run {s} -w {d} {rest} {p} {j} {c}'.format(t=threads, s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subflow.nf']))), d=workdir, rest=' '.join(argslist), p=' '.join("--{!s} {!s}".format(key,val) for (key,val) in params.items()), j=toolparams, c = '--CONDITION '+str.join(os.sep,condition))

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
            output, outerr = job.communicate()
            output = str.join('',output).rstrip()
            outerr = str.join('',outerr).rstrip()
            if output == '' and outerr == '' and job.poll() is not None:
                break
            if output and output != '':
                if any(x in output for x in ['ERROR','Error','error','Exception']) and not 'Workflow finished' in output or 'Execution complete' in output:
                    log.error(logid+'STOPPING: '+str(output))
                    log.info('PLEASE CHECK LOG AT LOGS/RunNextflow.log')
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid+str(output))
            if outerr and outerr != '':
                if not 'Workflow finished' in outerr and not 'Execution complete' in outerr and any(x in outerr for x in ['ERROR','Error','error','Exception']):
                    log.error(logid+'STOPPING: '+str(outerr))
                    log.info('PLEASE CHECK LOG AT LOGS/RunNextflow.log')
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid+str(outerr))
            if job.poll() is not None:
                break

        if job.returncode == 0:
            output, outerr = job.communicate()
            output = str.join('',output).rstrip()
            if output and output != '':
                log.info(logid+'JOB FINISHED: '+output)
            return job.poll()
        else:
            output, outerr = job.communicate()
            output = str.join('',output).rstrip()
            outerr = str.join('',outerr).rstrip()
            if outerr and outerr != '' or output and output != '':
                log.error(logid+'ERROR: '+outerr+output)
                log.info('PLEASE CHECK LOG AT LOGS/RunNextflow.log')
            job.kill()
            sys.exit('ERROR SIGNAL: '+str(job.returncode))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))
        sys.exit()

####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        knownargs=args[0]
        optionalargs=args[1:]

        log.setLevel(knownargs.loglevel)

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')
        log.debug(logid+str(log.handlers))

        run_nextflow(knownargs.configfile, knownargs.directory, knownargs.procs, knownargs.skeleton ,knownargs.loglevel, knownargs.clean, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# RunNextflow.py ends here
