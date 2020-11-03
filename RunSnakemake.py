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
# Last-Updated: Tue Nov  3 11:39:02 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 1101
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

import os, sys, json, shutil
import traceback as tb
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re

min_version("5.8.2")
scriptname=os.path.basename(__file__).replace('.py', '')

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
        config = load_configfile(configfile)
        if skeleton:
            for subdir in ['SubSnakes', 'RAW', 'FASTQ', 'LOGS', 'TMP']:
                makeoutdir(subdir)
            sys.exit('Skeleton directories created, please add files and rerun without --skeleton option')
        else:
            for subdir in ['SubSnakes', 'LOGS', 'TMP']:
                makeoutdir(subdir)

        subdir = 'SubSnakes'
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
            if '--profile' in optionalargs and 'nextsnakes/slurm' in optionalargs:
                makeoutdir('LOGS/SLURM')

        threads = min(int(config['MAXTHREADS']), procs) if 'MAXTHREADS' in config else procs

        if unlock:
            log.info(logid+'Unlocking directory')
            jobtorun = 'snakemake --unlock -j {t} -s {s} --configfile {c}'.format(t=threads, s=os.path.abspath(os.path.join('nextsnakes','workflows','header.smk')), c=configfile)
            log.info(logid+'UNLOCKING '+str(jobtorun))
            job = runjob(jobtorun)
            log.debug(logid+'JOB CODE '+str(job))

        preprocess = subworkflows = postprocess = None

        #Define workflow stages
        pre = ['QC','SRA']
        sub = ['QC','TRIMMING','MAPPING','DEDUP']
        post = ['COUNTING','UCSC','PEAKS','DE','DEU','DAS','ANNOTATE']

        wfs = config['WORKFLOWS'].split(',')

        if 'WORKFLOWS' in config:
            subworkflows = [x for x in wfs if x in sub]
            if len(subworkflows) == 0 or subworkflows[0] == '':
                subworkflows = None
            preprocess = [x for x in wfs if x in pre]
            if len(preprocess) == 0 or preprocess[0] == '':
                preprocess = None
            if subworkflows and 'MAPPING' in subworkflows and preprocess and 'QC' in preprocess:
                preprocess.remove('QC')
            postprocess = [x for x in wfs if x in post]
            if len(postprocess) == 0 or postprocess[0] == '':
                postprocess = None
        else:
            log.error('NO WORKFLOWS DEFINED, NOTHING TO DO!')
            sys.exit()

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

        if preprocess and 'SRA' in preprocess:
            if 'SRA' not in config:
                log.error(logid+'No configuration with key \'SRA\' for file download found. Nothing to do!')
            makeoutdir('FASTQ')
            makeoutdir('TMP')
            preprocess.remove('SRA')
            SAMPLES = download_samples(config)
            log.info(logid+'PRESAMPLES: '+str(SAMPLES))
            conditions = get_conditions(SAMPLES,config) #[x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
            log.info(logid+'PRECONDITIONS: '+str(conditions))
            for condition in conditions:
                subconf = NestedDefaultDict()
                subwork = 'SRA'
                listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                if listoftools is None:
                    log.warning(logid+'No entry fits condition '+str(condition)+' for preprocessing step '+str(subwork))
                    continue
                toolenv, toolbin = map(str,listoftools[0])
                if toolenv is None or toolbin is None:
                    continue
                subconf.update(listofconfigs[0])
                subname = toolenv+'.smk'
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.smk'))
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

                smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
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

        rawqc  = 'expand("QC/Multi/RAW/{condition}/multiqc_report.html", condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'

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
                        toolenv, toolbin = map(str,listoftools[i])
                        if toolenv is None or toolbin is None:
                            continue
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.smk'
                        log.debug(logid+'PREPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre_'+subwork,toolbin,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                            smkout.write('\n\n')

                        if subwork == 'QC':
                            subname = toolenv+'_raw.smk'

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write('rule themall:\n\tinput:\t'+rawqc+'\n\n')
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.smk'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())

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

        '''
        END OF PREPROCESSING, START OF PROCESSING
        '''

        if subworkflows:
            allmap = 'expand("MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))' if not 'DEDUP' in subworkflows else 'expand("MAPPED/{file}_mapped_sorted_unique_dedup.bam", file=samplecond(SAMPLES,config))'
            allqc  = 'expand("QC/Multi/{condition}/multiqc_report.html", condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'
            allrawqc  = 'expand("QC/Multi/RAW/{condition}/multiqc_report.html", condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'
            alltrimqc = 'expand("QC/Multi/TRIMMED_RAW/{condition}/multiqc_report.html",condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'
            alltrim = 'rule themall:\n    input: expand("TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz", file=samplecond(SAMPLES,config), read=["R1","R2"]) if paired == \'paired\' else expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=samplecond(SAMPLES,config))'
            alldedupqc = 'expand("QC/Multi/DEDUP_RAW/{condition}/multiqc_report.html",condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'
            alldedup = 'rule themall:\n    input: expand("DEDUP_FASTQ/{file}_{read}_dedup.fastq.gz", file=samplecond(SAMPLES,config), read=["R1","R2"]) if paired == \'paired\' else expand("DEDUP_FASTQ/{file}_dedup.fastq.gz", file=samplecond(SAMPLES,config))'
            alltrimdedupqc = 'expand("QC/Multi/DEDUP_TRIMMED_RAW/{condition}/multiqc_report.html", condition=str.join(os.sep,conditiononly(SAMPLES[0],config)))'

            log.info(logid+'STARTING PROCESSING')
            for condition in conditions:
                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.smk'))
                smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk'])))
                if os.path.exists(smko):
                    os.rename(smko,smko+'.bak')
                with open(smko, 'a') as smkout:
                    with open(smkf,'r') as smk:
                        for line in smk.readlines():
                            line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                            line = re.sub(condapath,'conda:  "../',line)
                            smkout.write(line)
                    smkout.write('\n\n')

                if 'QC' in subworkflows and 'QC' in config:
                    makeoutdir('QC')
                    if 'MAPPING' in subworkflows:
                        with open(smko, 'a') as smkout:
                            smkout.write('rule themall:\n\tinput:\t'+allmap+',\n\t\t'+allqc+'\n\n')
                    else:
                        if 'TRIMMING' in subworkflows and 'DEDUP' not in subworkflows:
                            with open(smko, 'a') as smkout:
                                smkout.write('rule themall:\n\tinput:\t'+alltrimqc+'\n\n')
                        elif 'TRIMMING' in subworkflows and 'DEDUP' in subworkflows:
                            with open(smko, 'a') as smkout:
                                smkout.write('rule themall:\n\tinput:\t'+alltrimdedupqc+'\n\n')
                        elif 'DEDUP' in subworkflows and 'TRIMMING' not in subworkflows:
                            with open(smko, 'a') as smkout:
                                smkout.write('rule themall:\n\tinput:\t'+alldedupqc+'\n\n')
                        else:
                            with open(smko, 'a') as smkout:
                                smkout.write('rule themall:\n\tinput:\t'+allrawqc+'\n\n')

                if 'MAPPING' in subworkflows and 'QC' not in subworkflows:
                    log.info(logid+'Mapping without QC for condition '+str(condition)+'!')
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write('rule themall:\n\tinput:\t'+allmap+'\n\n')
                        smkout.write('\n\n')

                if 'MAPPING' in subworkflows and 'TRIMMING' not in subworkflows:
                    log.info(logid+'Simulating read trimming as trimming is not part of the workflow for condition '+str(condition)+'!')
                    makeoutdir('TRIMMED_FASTQ')
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','simulatetrim.smk'))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')

                if 'TRIMMING' in subworkflows and 'QC' not in subworkflows and 'MAPPING' not in subworkflows:
                    log.info(logid+'Trimming without QC for condition'+str(condition)+'!')
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(alltrim+'\n')
                        smkout.write('\n\n')

                if 'DEDUP' in subworkflows and 'QC' not in subworkflows and 'TRIMMING' not in subworkflows and 'MAPPING' not in subworkflows:
                    log.info(logid+'DEDUP without QC for condition'+str(condition)+'!')
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(alldedup+'\n')
                        smkout.write('\n\n')

                subconf = NestedDefaultDict()
                for subwork in subworkflows:
                    log.debug(logid+'PREPARING '+str(subwork)+' '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        if toolenv is None or toolbin is None:
                            continue
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.smk'
                        log.debug(logid+'SUBWORKFLOW: '+str([subwork,toolenv,subname,condition, subsamples, subconf]))

                        if subwork == 'QC' and 'TRIMMING' in subworkflows and not 'MAPPING' in subworkflows:
                            if 'DEDUP' in subworkflows:
                                subname = toolenv+'_dedup_trim.smk'
                            else:
                                subname = toolenv+'_trim.smk'

                        if subwork == 'QC' and not 'TRIMMING' in subworkflows and not 'MAPPING' in subworkflows:
                            if 'DEDUP' in subworkflows:
                                subname = toolenv+'_dedup.smk'
                            else:
                                subname = toolenv+'_raw.smk'

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                if 'MAPPING' in subworkflows:
                    smkf = os.path.abspath(os.path.join('nextsnakes','workflows','mapping.smk'))
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')
                    if 'QC' in subworkflows:
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','multiqc.smk'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.smk'))
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
                jobtorun = 'snakemake -j {t} -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                log.info(logid+'RUNNING WORKFLOW '+str(jobtorun))
                job = runjob(jobtorun)
                log.debug(logid+'JOB CODE '+str(job))

        else:
            log.warning(logid+'No Workflows defined! Nothing to do, continuing with postprocessing!')


        if postprocess:
            for condition in conditions:
                for subwork in postprocess:
                    subconf = NestedDefaultDict()
                    if any(subwork == x for x in ['PEAKS', 'DE', 'DEU', 'DAS', 'DTU']):
                        continue

                    SAMPLES = get_samples_postprocess(config, subwork)
                    log.info(logid+'STARTING POSTPROCESSING '+str(subwork)+' WITH SAMPLES '+str(SAMPLES))
                    log.debug(logid+'POSTPROCESS: '+str(subwork)+' CONDITION: '+str(condition))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    log.debug(logid+str([listoftools,listofconfigs]))
                    if listoftools is None:
                        log.warning(logid+'No entry fits condition '+str(condition)+' for postprocessing step '+str(subwork))
                        continue

                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        if toolenv is None or toolbin is None:
                            continue
                        subconf.update(listofconfigs[i])
                        subname = toolenv+'.smk'
                        subsamples = list(set(sampleslong(subconf)))
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),subwork,toolbin,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix, 'loglevel=\''+loglevel+'\'', line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
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
            for analysis in ['PEAKS', 'DE', 'DEU', 'DAS', 'DTU']:
                if analysis in config and analysis in postprocess:

                    subwork = analysis
                    SAMPLES = get_samples_postprocess(config, subwork)
                    log.info(logid+'STARTING '+analysis+' Analysis '+' WITH SAMPLES '+str(SAMPLES))

                    subconf = NestedDefaultDict()
                    log.debug(logid+'SUBWORK: '+str(subwork)+' CONDITION: '+str(conditions))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, conditions)

                    if listoftools is None:
                        log.error(logid+'No entry fits condition '+str(conditions)+' for postprocessing step '+str(subwork))

                    for key in config[subwork]['TOOLS']:
                        log.info(logid+'... with Tool: '+key)
                        toolenv = key
                        toolbin = config[subwork]['TOOLS'][key]
                        subconf = NestedDefaultDict()
                        for i in listofconfigs:
                            if i is None:
                                continue
                            i[subwork+'ENV'] = toolenv
                            i[subwork+'BIN'] = toolbin
                        for i in range(len(listoftools)):
                            if listofconfigs[i] is None:
                                continue
                            subconf = merge_dicts(subconf,listofconfigs[i])

                        subname = toolenv+'.smk' if toolenv != 'edger' else toolenv+'_'+subwork+'.smk'
                        subsamples = sampleslong(subconf)
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','header.smk'))
                        smko = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                for line in smk.readlines():
                                    line = re.sub(logfix,'loglevel="'+loglevel+'"',line)
                                    line = re.sub(condapath,'conda:  "../',line)
                                    smkout.write(line)
                            smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows',subname))
                        with open(os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolenv,'subsnake.smk']))), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n')

                        smkf = os.path.abspath(os.path.join('nextsnakes','workflows','footer.smk'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())

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

    except Exception:
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
                if any(x in output for x in ['ERROR', 'Error', 'error', 'Exception']) and not 'Workflow finished' in output or 'Workflow finished' in outerr :
                    if outerr:
                        log.error(logid+'STOPPING: '+str(output)+'\n'+str(outerr))
                    else:
                        log.error(logid+'STOPPING: '+str(output))
                    log.info('PLEASE CHECK LOG AT LOGS/RunSnakemake.log')
                    job.kill()
                    sys.exit()
                else:
                    log.info(logid+str(output))
            if outerr and outerr != '':
                if not 'Workflow finished' in outerr and not 'Nothing to be done' in outerr and not 'Workflow finished' in output and any(x in outerr for x in ['ERROR', 'Error', 'error', 'Exception']):
                    log.error(logid+'STOPPING: '+str(outerr))
                    log.info('PLEASE CHECK LOG AT LOGS/RunSnakemake.log')
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
                log.info('PLEASE CHECK LOG AT LOGS/RunSnakemake.log')
            job.kill()
            sys.exit('ERROR SIGNAL: '+str(job.returncode))

    except Exception:
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

        run_snakemake(knownargs.configfile, knownargs.debug_dag, knownargs.filegraph, knownargs.directory, knownargs.use_conda, knownargs.procs, knownargs.skeleton, knownargs.loglevel, knownargs.unlock, optionalargs[0])

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# RunSnakemake.py ends here
