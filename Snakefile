import glob, os, sys, inspect, snakemake, json, shutil
from collections import defaultdict
import traceback as tb
import subprocess
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
    postprocess = config['POSTPROCESSING'].split(',')  # we keep this separate because not all postprocessing steps need extra configuration

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

    log.debug(logid+'SAMPLES: '+str(SAMPLES))
    conditions = [x.split(os.sep) for x in list(set([os.path.dirname(x) for x in samplecond(SAMPLES,config)]))]
    log.debug(logid+'CONDITIONS: '+str(conditions))

    for condition in conditions:
        if 'MAPPING' in subworkflows:
            smkf = os.path.abspath(os.path.join('snakes','workflows','header.smk'))
            with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                with open(smkf,'r') as smk:
                    smkout.write(smk.read())

            smkf = os.path.abspath(os.path.join('snakes','workflows','mapping.smk'))
            with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                with open(smkf,'r') as smk:
                    smkout.write(smk.read())

        if 'QC' in subworkflows:
            smkf = os.path.abspath(os.path.join('snakes','workflows','multiqc.smk'))
            with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                with open(smkf,'r') as smk:
                    smkout.write(smk.read())

        subconf = NestedDefaultDict()
        for subwork in subworkflows:
            listoftools, listofconfigs = create_subworkflow(config, subwork, [condition])
            for i in range(0,len(listoftools)):
                toolenv, toolbin = map(str,listoftools[i])
                subconf.update(listofconfigs[i])
                subsamples = list(set(sampleslong(subconf)))
                paired = ''
                if checkpaired(SAMPLES, config):
                    paired = '_paired'
                subname = toolenv+paired+'.smk'
                log.debug(logid+'SUBWORKFLOW: '+str([toolenv,subname,condition, subsamples, subconf]))

                smkf = os.path.abspath(os.path.join('snakes','workflows',subname))
                with open('_'.join(['_'.join(condition),'subsnake.smk']), 'a') as smkout:
                    with open(smkf,'r') as smk:
                        smkout.write(smk.read())

        with open('_'.join(['_'.join(condition),'subconfig.json']), 'a') as confout:
            json.dump(subconf, confout)

    for condition in conditions:
        jobtorun = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory ../. --printshellcmds --debug-dag'.format(t=MAXTHREAD,s='_'.join(['_'.join(condition),'subsnake.smk']),c='_'.join(['_'.join(condition),'subconfig.json']))
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

                jobtorun  = 'snakemake -j {t} --use-conda -s {s} --configfile {c} --directory . --printshellcmds --debug-dag'.format(t=MAXTHREAD,s='_'.join(['_'.join(condition),toolbin,'subsnake.smk']),c='_'.join(['_'.join(condition),toolbin,'subconfig.json']))
                o = subprocess.run(jobtorun, shell=True, check=True, universal_newlines=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                if o.stderr:
                    log.error(o.sterr)

except Exception as err:
    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))
