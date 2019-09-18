import glob, os, sys, inspect, snakemake, json, shutil
from collections import defaultdict
import traceback as tb
from snakemake.utils import validate, min_version
min_version("5.5.2")

###snakemake -n -j 20 --use-conda -s Workflow/Snakefile
###--configfile Workflow/configs/config_example.json --directory ${PWD}
###--printshellcmds &> LOGS/run.log

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from lib.Collection import *
from lib.Logger import *

try:
    log = setup_logger(name='Snakemake', log_file='LOGS/Snakemake.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')

    subworkflows = config['WORKFLOWS'].split(',')

    try:
        all([config[x] for x in subworkflows])
    except KeyError:
        log.error('Not all subworkflows have configuration in the config file')

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
    except KeyError:
        CLIP=''

    conditions = list(set([join(os.sep,samplecond(x,config)) for x in SAMPLES]))

    for subwork in subworkflows:
        listoftools, listofconfigs = create_subworkflow(config, subwork, conditions)
        for i in range(0,len(listoftools)):
            toolenv, toolbin = listoftools[i]
            subconf = listofconfigs[i]
            with open('_'.join(['subconfig',toolenv,'_'.join(conditions),'subworkflow.json']), 'w') as outfile:
                json.dump(tempconf, outfile)
            subworkflow sampleqc:
                snakefile: 'workflows/'+toolenv+'.smk'
                configfile: '_'.join([toolenv,'_'.join(conditions),'subworkflow.json'])
                workdir: '.'

        rule all:
            input: sampleqc(expand("DONE/{file}"+toolenv, file=samplecond(SAMPLES,config)))

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))
