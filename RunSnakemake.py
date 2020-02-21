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
# Last-Updated: Fri Feb 21 11:09:21 2020 (+0100)
#           By: Joerg Fallmann
#     Update #: 632
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
min_version("5.8.2")

from lib.Collection import *
from lib.Logger import *        # Switch to snakemake logger?
scriptname=os.path.basename(__file__)

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

def run_snakemake (configfile, debugdag, filegraph, workdir, useconda, procs, skeleton, unlock=None, optionalargs=None):
    try:
        logid = scriptname+'.run_snakemake: '
        if skeleton:
            for subdir in ['SubSnakes', 'RAW', 'GENOMES', 'FASTQ', 'LOGS']:  # Add RAW for nanopore preprocessing
                makeoutdir(subdir)
            sys.exit('Skeleton directories created, please add files and rerun without --skeleton option')
        else:
            for subdir in ['SubSnakes', 'LOGS']:  # Add RAW for nanopore preprocessing
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

        if unlock:
            log.info(logid+'Unlocking directory')
            jobtorun = 'snakemake --unlock -s {s} --configfile {c}'.format(s=os.path.abspath(os.path.join('snakes','workflows','header.smk')), c=configfile)
            log.info(logid+'RUNNING '+str(jobtorun))
            o = runjob(jobtorun)
            if o.stdout:
                log.info(o.stdout)
                if any(x in o.stdout for x in ['ERROR','Error','error']):
                    sys.exit(o.stdout)

            if o.stderr:
                log.error(o.stderr)
                if any(x in o.stderr for x in ['ERROR','Error','error']):
                    sys.exit(o.stderr)

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

        log.debug(logid+str([preprocess,subworkflows,postprocess]))

        threads = min(int(config['MAXTHREADS']), procs) if 'MAXTHREADS' in config else procs

        # CLEANUP
        #oldsmk = os.path.abspath(os.path.join(subdir,'*_subsnake.smk'))
        #oldcnf = os.path.abspath(os.path.join(subdir,'*_subconfig.json'))
        #for oldfile in glob.glob(oldsmk):
        #    os.rename(oldfile,oldfile+'.bak')
        #    log.warning(logid+'Found old snakemake file'+oldfile+', was moved to '+oldfile+'.bak')
        #for oldfile in glob.glob(oldcnf):
        #    os.rename(oldfile,oldfile+'.bak')
        #    log.warning(logid+'Found old config file'+oldfile+', was moved to '+oldfile+'.bak')

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

        #subworkflows.extend(postprocess)  # Concatenate to get the full list of steps to process
        log.debug(logid+'WORKFLOWS: '+str(subworkflows))

        SAMPLES=[os.path.join(x) for x in list(set(samples(config)))]
        check = [os.path.join('FASTQ',x+'*.fastq.gz') for x in SAMPLES]
        log.debug(logid+'SAMPLECHECK: '+str(check))
        sampletest = [glob.glob(os.path.abspath(x)) for x in check][0]
        log.debug(logid+'SAMPLETEST: '+str(sampletest))
        if len(sampletest) < 1:
            SAMPLES = [os.path.join(x) for x in sampleslong(config)]
            log.debug(logid+'SAMPLES_LONG: '+str(SAMPLES))
            check = [os.path.join('FASTQ',str(x)+'*.fastq.gz') for x in SAMPLES]
            log.debug(logid+'SAMPLECHECK_LONG: '+str(check))
            log.debug(logid+'SAMPLECHECK_LONG: '+str([glob.glob(os.path.abspath(x)) for x in check][0]))
            sampletest = [glob.glob(os.path.abspath(x)) for x in check][0]
            log.debug(logid+'SAMPLETEST_LONG: '+str(sampletest))
            if len(sampletest) < 1:
                log.error(logid+'No samples found, please check config file')
                sys.exit()

        log.info(logid+'Working on SAMPLES: '+str(SAMPLES))
        conditions = [x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
        log.info(logid+'CONDITIONS: '+str(conditions))

        condapath=re.compile(r'conda:\s+"')

        if preprocess:
            log.info(logid+'STARTING PREPROCESSING')
            if 'QC' in preprocess and 'QC' in config:
                makeoutdir('QC')
            for condition in conditions:
                subconf = NestedDefaultDict()
                for subwork in preprocess:

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
                        smko = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre'+subwork,toolbin,'subsnake.smk'])))
                        if os.path.exists(smko):
                            os.rename(smko,smko+'.bak')
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                                smkout.write('\n\n')

                        if subwork == 'QC':
                            subname = toolenv+'_raw.smk'
                        #    with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre'+subwork,toolbin,'subsnake.smk']))), 'a') as smkout:
                        #        smkf = os.path.abspath(os.path.join('snakes','workflows','premultiqc.smk'))
                        #        with open(smkf,'r') as smk:
                        #            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        #            smkout.write('\n\n')

                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                                smkout.write('\n\n')

                        confo = os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre'+subwork,toolbin,'subconfig.json'])))
                        if os.path.exists(confo):
                            os.rename(confo,confo+'.bak')
                        with open(confo, 'a') as confout:
                            json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre'+subwork,toolbin,'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'pre'+subwork,toolbin,'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                        log.info(logid+'RUNNING '+str(jobtorun))
                        o = runjob(jobtorun)
                        if o.stdout:
                            log.info(o.stdout)
                            if not 'Workflow finished, no error' in o.stdout or 'Exception' in o.stdout:
                                sys.exit(o.stdout)

                        if o.stderr:
                            log.error(o.stderr)
                            if any(x in o.stderr for x in ['ERROR','Error','error','Exception']):
                                sys.exit(o.stderr)
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
                        smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                    smkout.write('\n\n')

                if 'QC' in subworkflows and 'QC' in config:
                    makeoutdir('QC')
                    if 'MAPPING' in subworkflows:
                        with open(smko, 'a') as smkout:
                            smkout.write('rule all:\n\tinput: expand("DONE/{file}_mapped",file=samplecond(SAMPLES,config))\n\n')

                        smkf = os.path.abspath(os.path.join('snakes','workflows','multiqc.smk'))
                        with open(smko, 'a') as smkout:
                            with open(smkf,'r') as smk:
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
                with open(confout, 'a') as confout:
                    json.dump(subconf, confout)

            for condition in conditions:
                log.info(logid+'Starting workflows for condition '+str(condition))
                jobtorun = 'snakemake -j {t} -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))),d=workdir,rest=' '.join(argslist))
                log.info(logid+'RUNNING WORKFLOW '+str(jobtorun))
                o = runjob(jobtorun)
                if o.stdout:
                    log.info(o.stdout)
                    if not 'Workflow finished, no error' in o.stdout or 'Exception' in o.stdout:
                        sys.exit(o.stdout)

                if o.stderr:
                    log.error(o.stderr)
                    if any(x in o.stderr for x in ['ERROR','Error','error','Exception']):
                        sys.exit(o.stderr)
        else:
            log.warning(logid+'No subworkflows defined! Nothing to do!')

        if postprocess:
            log.info(logid+'STARTING POSTPROCESSING')

            if 'PEAKS' in config and 'PEAKS' in postprocess:
                CLIP = checkclip(SAMPLES, config)
                log.info(logid+'Running Peak finding for '+CLIP+' protocol')

            for condition in conditions:
                subconf = NestedDefaultDict()
                for subwork in postprocess:
                    if subwork == 'DE' or subwork == 'DEU' or subwork == 'DAS':
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
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
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
                        o = runjob(jobtorun)
                        if o.stdout:
                            log.info(o.stdout)
                            if not 'Workflow finished, no error' in o.stdout or 'Exception' in o.stdout:
                                #if any(x in o.stdout for x in ['ERROR','Error','error']):
                                sys.exit(o.stdout)

                        if o.stderr:
                            log.error(o.stderr)
                            if any(x in o.stderr for x in ['ERROR','Error','error','Exception']):
                                sys.exit(o.stderr)

            #THIS SECTION IS FOR DE, DEU, DAS ANALYSIS, WE USE THE CONDITIONS TO MAKE PAIRWISE COMPARISONS
            for analysis in ['DE', 'DEU', 'DAS']:
                if analysis in config and analysis in postprocess:
                    log.info(logid+'STARTING '+analysis+' Analysis')
                    subwork = analysis
                    subconf = NestedDefaultDict()
                    log.debug(logid+'SUBWORK: '+str(subwork)+' CONDITION: '+str(conditions))

                    listoftools, listofconfigs = create_subworkflow(config, subwork, conditions)
                    listoftoolscount, listofconfigscount = create_subworkflow(config, 'COUNTING', conditions)
                    if listoftools is None or listoftoolscount is None:
                        log.error(logid+'No entry fits condition '+str(conditions)+' for postprocessing step '+str(subwork)+' or COUNTING not configured')

                    toolenv, toolbin = map(str,listoftools[0])
                    countenv, countbin = map(str,listoftoolscount[0])

                    subconf = listofconfigs[0]
                    for x in range(1,len(listofconfigs)):
                        subconf = merge_dicts(subconf,listofconfigs[x])
                    for x in range(0,len(listofconfigscount)):
                        subconf = merge_dicts(subconf,listofconfigscount[x])

                    subname = toolenv+'.smk'
                    subsamples = sampleslong(subconf)
                    log.debug(logid+'POSTPROCESS: '+str([toolenv,subname, subsamples, subconf]))
                    smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                    smko = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolbin,'subsnake.smk'])))
                    if os.path.exists(smko):
                        os.rename(smko,smko+'.bak')
                    with open(smko, 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')
                    smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                    with open(os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolbin,'subsnake.smk']))), 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                    confo = os.path.abspath(os.path.join(subdir,'_'.join([subwork,toolbin,'subconfig.json'])))
                    if os.path.exists(confo):
                        os.rename(confo,confo+'.bak')
                    with open(confo, 'a') as confout:
                        json.dump(subconf, confout)

                    jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join([subwork,toolbin,'subsnake.smk'])]))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join([subwork,toolbin,'subconfig.json'])]))),d=workdir,rest=' '.join(argslist))
                    log.info(logid+'RUNNING '+str(jobtorun))
                    o = runjob(jobtorun)
                    if o.stdout:
                        log.info(o.stdout)
                        if not 'Workflow finished, no error' in o.stdout or 'Exception' in o.stdout:
                            #if any(x in o.stdout for x in ['ERROR','Error','error']):
                            sys.exit(o.stdout)

                    if o.stderr:
                        log.error(o.stderr)
                        if any(x in o.stderr for x in ['ERROR','Error','error','Exception']):
                            sys.exit(o.stderr)

        else:
            log.warning(logid+'No postprocessing steps defined! Nothing to do!')

        log.info('Workflows executed without error!')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(''.join(tbe.format()))


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
        log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')

        run_snakemake(knownargs.configfile, knownargs.debug_dag, knownargs.filegraph, knownargs.directory, knownargs.use_conda, knownargs.procs, knownargs.skeleton, knownargs.unlock, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

# RunSnakemake.py ends here
