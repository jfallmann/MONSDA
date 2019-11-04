#!/usr/bin/env python3

import glob, os, sys, inspect, json, shutil
from collections import defaultdict
import traceback as tb
#import snakemake
from snakemake import load_configfile
from snakemake.utils import validate, min_version
import argparse
import subprocess
import re
min_version("5.7.1")

#cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
#if cmd_subfolder not in sys.path:
#    sys.path.insert(0, cmd_subfolder)

from lib.Collection import *
from lib.Logger import *
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
    parser.add_argument("-v", "--loglevel", type=str, default='INFO', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_known_args()

def run_snakemake (configfile, debugdag, filegraph, workdir, useconda, procs, unlock=None, optionalargs=None):
    try:
        for subdir in ['SubSnakes', 'GENOMES', 'FASTQ', 'TRIMMED_FASTQ', 'QC', 'LOGS']:
            makeoutdir(subdir)

        subdir = 'SubSnakes'
        config = load_configfile(configfile)
        argslist = list()
        if useconda:
            argslist.append("--use-conda")
        if debugdag:
            argslist.append("--debug-dag")
        if filegraph:
            argslist.append("--filegraph|dot|display")
        if optionalargs and len(optionalargs) > 0:
            argslist.extend(optionalargs)

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

        subworkflows = config['WORKFLOWS'].split(',')
        if len(subworkflows) == 0 or subworkflows[0] == '':
            subworkflows = None

        postprocess = config['POSTPROCESSING'].split(',') # we keep this separate because not all postprocessing steps need extra configuration
        if len(postprocess) == 0 or postprocess[0] == '':
            postprocess = None

        threads = min(int(config['MAXTHREADS']), procs)

        # CLEANUP
        oldsmk = os.path.abspath(os.path.join(subdir,'*_subsnake.smk'))
        oldcnf = os.path.abspath(os.path.join(subdir,'*_subconfig.json'))
        for oldfile in glob.glob(oldsmk):
            os.rename(oldfile,oldfile+'.bak')
            log.warning(logid+'Found old snakemake file'+oldfile+', was moved to '+oldfile+'.bak')
        for oldfile in glob.glob(oldcnf):
            os.rename(oldfile,oldfile+'.bak')
            log.warning(logid+'Found old config file'+oldfile+', was moved to '+oldfile+'.bak')

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

        SAMPLES=list(set(samples(config)))
        if os.path.exists(SAMPLES[0]) is False:
            SAMPLES=list(set(sampleslong(config)))

        log.info(logid+'Working on SAMPLES: '+str(SAMPLES))
        conditions = [x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
        log.info(logid+'CONDITIONS: '+str(conditions))

        condapath=re.compile(r'conda:\s+"')

        if subworkflows:
            for condition in conditions:
                smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                    smkout.write('\n\n')

                if 'QC' in subworkflows and config['QC']['RUN'] == "ON":
                    if 'MAPPING' in subworkflows:
                        with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                            smkout.write('rule all:\n\tinput: expand("DONE/{file}_mapped",file=samplecond(SAMPLES,config))\n\n')

                    smkf = os.path.abspath(os.path.join('snakes','workflows','multiqc.smk'))
                    with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')

                if 'MAPPING' in subworkflows and ('TRIMMING' not in subworkflows or ('TRIMMING' in config and config['TRIMMING']['RUN'] == "OFF")):
                    log.info(logid+'Simulating read trimming as trimming was set to OFF or is not part of the workflow!')
                    smkf = os.path.abspath(os.path.join('snakes','workflows','simulatetrim.smk'))
                    with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                        smkout.write('\n\n')

                subconf = NestedDefaultDict()
                for subwork in subworkflows:
                    if 'RUN' in config[subwork] and config[subwork]['RUN'] == "OFF":
                        log.info(logid+'Workflowstep '+subwork+' set to OFF, will be skipped')
                        continue
                    else:
                        listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                        for i in range(0,len(listoftools)):
                            toolenv, toolbin = map(str,listoftools[i])
                            subconf.update(listofconfigs[i])
                            subsamples = list(set(sampleslong(subconf)))
                            subname = toolenv+'.smk'
                            log.debug(logid+'SUBWORKFLOW: '+str([subwork,toolenv,subname,condition, subsamples, subconf]))

                            smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                            with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                                with open(smkf,'r') as smk:
                                    smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                                smkout.write('\n\n')

                if 'MAPPING' in subworkflows:
                    smkf = os.path.abspath(os.path.join('snakes','workflows','mapping.smk'))
                    with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))), 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                            smkout.write('\n\n')

                with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))), 'a') as confout:
                    json.dump(subconf, confout)

            for condition in conditions:
                log.info(logid+'Starting runs for condition '+str(condition))
                jobtorun = 'snakemake -j {t} -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),'subconfig.json']))),d=workdir,rest=' '.join(argslist))
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
            log.warning(logid+'No subworkflows defined! Nothing to do!')

        if postprocess:
            log.info(logid+'Starting runs for postprocessing')

            if 'PEAKS' in config and 'PEAKS' in postprocess:
                CLIP = checkclip(SAMPLES, config)
                log.info(logid+'Running Peak finding for '+CLIP+' protocol')

            for condition in conditions:
                subconf = NestedDefaultDict()
                for subwork in postprocess:
                    log.debug(logid+'SUBWORK: '+str(subwork))
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    log.debug(logid+str([listoftools,listofconfigs]))
                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subname = toolenv+'.smk'
                        subsamples = list(set(sampleslong(subconf)))
                        log.debug(logid+'POSTPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))
                        smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
                        with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),toolbin,'subsnake.smk']))), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                                smkout.write('\n\n')
                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),toolbin,'subsnake.smk']))), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(re.sub(condapath,'conda:  "../',smk.read()))
                                smkout.write('\n\n')
                        listoftools, listofconfigs = create_subworkflow(config, "ANNOTATE", [condition])
                        subconf.update(listofconfigs[i])
                        with open(os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),toolbin,'subconfig.json']))), 'a') as confout:
                            json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds --show-failed-logs {rest}'.format(t=threads,s=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),toolbin,'subsnake.smk']))),c=os.path.abspath(os.path.join(subdir,'_'.join(['_'.join(condition),toolbin,'subconfig.json']))),d=workdir,rest=' '.join(argslist))
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
        log.addHandler(logging.StreamHandler(sys.stdout))  # streamlog

        MIN_PYTHON = (3,7)
        if sys.version_info < MIN_PYTHON:
            log.error("This script requires Python version >= 3.7")
            sys.exit("This script requires Python version >= 3.7")
        log.info(logid+'Running '+scriptname+' on '+str(knownargs.procs)+' cores')

        run_snakemake(knownargs.configfile, knownargs.debug_dag, knownargs.filegraph, knownargs.directory, knownargs.use_conda, knownargs.procs, knownargs.unlock, optionalargs[0])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
