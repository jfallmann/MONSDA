import glob
import os
import sys
import inspect
import snakemake
import json
import shutil
import logging
import tempfile
import traceback as tb
from collections import defaultdict
from itertools import combinations

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../nextsnakes/lib"),os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"nextsnakes/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from Collection import *
from Logger import *

loglevel="INFO"

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace('.py','')
    if any(x in scriptname for x in ['RunSnakemake','Configurator']):
        log = logging.getLogger(scriptname)
    else:
        log = logging.getLogger('snakemake')
        for handler in log.handlers[:]:
            handler.close()
            log.removeHandler(handler)
    handler = logging.FileHandler('LOGS/RunSnakemake.log', mode='a')
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s',datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s',datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    lvl = loglevel
    log.setLevel(lvl)

except Exception as err:
    log = setup_logger(name='RunSnakemake.header', log_file='LOGS/RunSnakemake.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level=loglevel, filemode='a')
    #log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=loglevel)

    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))

logid = 'header.smk: '

#Config Parsing
BINS = config["BINS"]
MAXTHREAD = int(config["MAXTHREADS"])
SAMPLES = [os.path.join(x) for x in sampleslong(config)]
SETTINGS = keysets_from_dict(config["SAMPLES"])[0]
SETS = os.sep.join(SETTINGS)
log.debug(logid+'SETS: '+SETS)

if len(SAMPLES) < 1:
    log.error(logid+'No samples found, please check config file')
    sys.exit(logid+'ERROR: No samples found, please check config file')

log.info(logid+'Working on SAMPLES: '+str(SAMPLES))

paired = checkpaired([SAMPLES[0]], config)
if paired == 'paired':
    log.info('RUNNING SNAKEMAKE IN PAIRED READ MODE')

stranded = checkstranded([SAMPLES[0]], config)
if stranded != '':
    log.info('RUNNING SNAKEMAKE WITH STRANDEDNESS '+str(stranded))


#MAPPING Variables
if 'MAPPING' in config:
    MAPCONF = subDict(config['MAPPING'], SETTINGS)
    log.debug(logid+'MAPPINGCONFIG: '+str(SETTINGS)+'\t'+str(MAPCONF))
    REFERENCE = MAPCONF['REFERENCE']
    REFDIR = str(os.path.dirname(REFERENCE))
    ANNO = MAPCONF.get('ANNOTATION')
    MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')
    log.debug(logid+'INDEX: '+str(MAPCONF['INDEX']))
    UIDX = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))
    INDICES = MAPCONF['INDEX'].split(',') if 'INDEX' in MAPCONF else list(UIDX)
    INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
    PREFIX = MAPCONF['PREFIX'] if 'PREFIX' in MAPCONF else ''

    if len(INDICES) > 1:
        if str(os.path.abspath(INDICES[1])) not in UIDX:
            INDEX2 = str(os.path.abspath(INDICES[1]))
        else:
            INDEX2 = str(os.path.abspath(INDICES[1]))+'_idx'
    else:
        INDEX2 = ''

    log.debug(logid+'REF: '+'\t'.join([REFERENCE,REFDIR,INDEX,str(INDEX2)]))

#Peak Calling Variables
if 'PEAKS' in config:
    PEAKCONF = subDict(config['PEAKS'], SETTINGS)
    REFERENCE = PEAKCONF['REFERENCE']
    REFDIR = str(os.path.dirname(REFERENCE))
    ANNOPEAK = PEAKCONF.get('ANNOTATION')
    CLIP = checkclip(SAMPLES, config)
    log.info(logid+'Running Peak finding for '+CLIP+' protocol')
    peakcallconf = tool_params(SAMPLES[0],None,config,'PEAKS')['OPTIONS'][0]
    try:
        all([x in peakcallconf for x in ['MINPEAKRATIO', 'PEAKDISTANCE', 'PEAKWIDTH', 'PEAKCUTOFF', 'MINPEAKHEIGHT', 'USRLIMIT']])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error('Not all required peak finding options defined in config!\n'+''.join(tbe.format()))

    MINPEAKRATIO = peakcallconf['MINPEAKRATIO']
    PEAKDISTANCE = peakcallconf['PEAKDISTANCE']
    PEAKWIDTH = peakcallconf['PEAKWIDTH']
    PEAKCUTOFF = peakcallconf['PEAKCUTOFF']
    MINPEAKHEIGHT = peakcallconf['MINPEAKHEIGHT']
    USRLIMIT = peakcallconf['USRLIMIT']

    if 'PREPROCESS' in peakcallconf:
        PREPROCESS = ' '.join("{!s} {!s}".format(key,val) for (key,val) in peakcallconf['PREPROCESS'].items())
    else:
        PREPROCESS = ''
