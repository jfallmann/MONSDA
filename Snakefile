import glob, os, sys, inspect, snakemake, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake.utils import validate, min_version
min_version("5.5.2")

###snakemake -n -j 20 --use-conda -s Workflow/Snakefile
###--configfile Workflow/configs/config_example.json --directory ${PWD}
###--printshellcmds &> LOGS/run.log

from lib.Logger import *
from lib.Collection import *

try:
    log = setup_logger(name='Snakemake', log_file='LOGS/Snakemake.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')
    logid = 'Main Snakefile: '
    subworkflows = config['WORKFLOWS'].split(',')
    #postprocess = config['POSTPROCESS'].split(',')  # we keep this separate because not all postprocessing steps need extra configuration

    try:
        all([config[x] for x in subworkflows])
    except KeyError:
        log.warning('Not all required subworkflows have configuration in the config file')

    #subworkflows.extend(postprocess)  # Concatenate to get the full list of steps to process
    log.debug(subworkflows)

    QC=config["QC"]
    REFERENCE=config["REFERENCE"]
    GENOME=config["GENOME"]
    NAME=config["NAME"]
    BINS=config["BINS"]
    MAXTHREAD=int(config["MAXTHREADS"])
    SOURCE=sources(config)
    SAMPLES=list(set(samples(config)))
    if os.path.exists(SAMPLES[0]) is False:
        SAMPLES=list(set(sampleslong(config)))
    try:
        CLIP=config["CLIP"]
    except:
        CLIP=''

    log.debug(logid+str(SAMPLES))
    conditions = [x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
    log.debug(logid+str(conditions))

    for subwork in subworkflows:
        listoftools, listofconfigs = create_subworkflow(config, subwork, conditions)
        for i in range(0,len(listoftools)):
            toolenv, toolbin = map(str,listoftools[i])
            subconf = listofconfigs[i]
            subname = toolenv+'.smk'
            log.debug([toolenv,subname,conditions])
            if not (os.path.islink(subname) or os.path.isfile(subname)):
                os.symlink(os.path.join('snakes','workflows',subname),subname)
            with open('_'.join(['subconfig',toolenv,'_'.join(conditions[i]),'subworkflow.json']), 'w') as outfile:
                json.dump(subconf, outfile)
            subworkflow do_work:
                workdir: '../.'
                snakefile: subname
                configfile: '_'.join(['subconfig',toolenv,'_'.join(conditions[i]),'subworkflow.json'])

            rule all:
                input: do_work(expand("DONE/{file}"+toolenv, file=samplecond(SAMPLES,config)))

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))
