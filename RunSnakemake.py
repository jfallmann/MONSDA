#!/usr/bin/env python3

import glob, os, sys, inspect, snakemake, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake.utils import validate, min_version
from snakemake import load_configfile
import argparse
import subprocess
min_version("5.5.2")

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
    parser.add_argument("-d", "--directory", type=str, default='', help='Directory to work in')
    parser.add_argument("-u", "--use-conda", action="store_true", help='Should conda be used')
    parser.add_argument("-j", "--procs", type=int, default=1, help='Number of parallel processed to start snakemake with, capped by MAXTHREADS in config!')
    parser.add_argument("-v", "--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

def run_snakemake (configfile, debugdag, workdir, useconda, procs):
    try:
        config = load_configfile(configfile)
        if useconda :
            useconda = "--use-conda"
        else:
            useconda = ''
        if debugdag :
            debugdag = "--debug-dag"
        else:
            debugdag = ''

        subworkflows = config['WORKFLOWS'].split(',')
        postprocess = config['POSTPROCESSING'].split(',')  # we keep this separate because not all postprocessing steps need extra configuration

        threads = min(int(config['MAXTHREADS']), procs)

        # CLEANUP
        oldsmk = os.path.abspath(os.path.join('*_subsnake.smk'))
        oldcnf = os.path.abspath(os.path.join('*_subconfig.json'))
        for oldfile in glob.glob(oldsmk):
            log.warning(logid+'Found old snakemake file'+oldfile+', will be moved to '+oldfile+'.bak')
            os.rename(oldfile,oldfile+'.bak')
        for oldfile in glob.glob(oldcnf):
            log.warning(logid+'Found old config file'+oldfile+', will be moved to '+oldfile+'.bak')
            os.rename(oldfile,oldfile+'.bak')

        try:
            all([config[x] for x in subworkflows])
        except KeyError:
            log.warning(logid+'Not all required subworkflows have configuration in the config file')

        #subworkflows.extend(postprocess)  # Concatenate to get the full list of steps to process
        log.debug(logid+'WORKFLOWS: '+str(subworkflows))

        SAMPLES=list(set(samples(config)))
        if os.path.exists(SAMPLES[0]) is False:
            SAMPLES=list(set(sampleslong(config)))
        try:
            CLIP=config["CLIP"]
        except:
            CLIP=''

        log.debug(logid+'SAMPLES: '+str(SAMPLES))
        conditions = [x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
        log.debug(logid+'CONDITIONS: '+str(conditions))

        for condition in conditions:
            smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
            with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                with open(smkf,'r') as smk:
                    smkout.write(smk.read())
                smkout.write('\n\n')

            if 'QC' in subworkflows:
                smkf = os.path.abspath(os.path.join('snakes','workflows','multiqc.smk'))
                with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())
                    smkout.write('\n\n')

            if 'MAPPING' in subworkflows:
                subconf = NestedDefaultDict()
                for subwork in subworkflows:
                    listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                    for i in range(0,len(listoftools)):
                        toolenv, toolbin = map(str,listoftools[i])
                        subconf.update(listofconfigs[i])
                        subsamples = list(set(sampleslong(subconf)))
                        subname = toolenv+'.smk'
                        log.debug(logid+'SUBWORKFLOW: '+str([toolenv,subname,condition, subsamples, subconf]))

                        smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                        with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                            with open(smkf,'r') as smk:
                                smkout.write(smk.read())
                            smkout.write('\n\n')

                smkf = os.path.abspath(os.path.join('snakes','workflows','mapping.smk'))
                with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())
                    smkout.write('\n\n')

            with open('_'.join(['_'.join(condition),'subconfig.json']), 'a') as confout:
                json.dump(subconf, confout)

        for condition in conditions:
            jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds {g}'.format(t=threads,s='_'.join(['_'.join(condition),'subsnake.smk']),c='_'.join(['_'.join(condition),'subconfig.json']),d=workdir,g=debugdag)
            o = subprocess.run(jobtorun, shell=True, check=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if o.stderr:
                log.error(o.stderr)

        for condition in conditions:
            subconf = NestedDefaultDict()
            for subwork in postprocess:
                listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
                for i in range(0,len(listoftools)):
                    toolenv, toolbin = map(str,listoftools[i])
                    subconf.update(listofconfigs[i])
                    subname = toolenv+'.smk'
                    subsamples = list(set(sampleslong(subconf)))
                    log.debug(logid+'POSTPROCESS: '+str([toolenv,subname,condition, subsamples, subconf]))

                    smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                    with open('_'.join(['_'.join(condition),toolbin,'subsnake.smk']), 'a') as smkout:
                        with open(smkf,'r') as smk:
                            smkout.write(smk.read())
                    with open('_'.join(['_'.join(condition),toolbin,'subconfig.json']), 'a') as confout:
                        json.dump(subconf, confout)

                        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory {d} --printshellcmds {g}'.format(t=MAXTHREAD,s='_'.join(['_'.join(condition),'subsnake.smk']),c='_'.join(['_'.join(condition),'subconfig.json']),d=workdir,g=debugdag)
                    o = subprocess.run(jobtorun, shell=True, check=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if o.stderr:
                        log.error(o.sterr)
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
        log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
        log.addHandler(logging.StreamHandler(sys.stdout))  # streamlog
        log.info(logid+'Running '+scriptname+' on '+str(args.procs)+' cores')

        run_snakemake(args.configfile, args.debug_dag, args.directory, args.use_conda, args.procs)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
