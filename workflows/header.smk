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
import re

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../NextSnakes/NextSnakes"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"NextSnakes/NextSnakes"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../NextSnakes"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"NextSnakes")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from NextSnakes.Logger import *
from NextSnakes.Params import *


loglevel="INFO"

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace('.py','')
    if any(x in scriptname for x in ['NextSnakes','Configurator']):
        log = logging.getLogger(scriptname)
    else:
        log = logging.getLogger('snakemake')
        for handler in log.handlers[:]:
            handler.close()
            log.removeHandler(handler)
    handler = logging.FileHandler('LOGS/NextSnakes.log', mode='a')
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    lvl = loglevel
    log.setLevel(lvl)

except Exception as err:
    log = setup_logger(name='NextSnakes.header', log_file='LOGS/NextSnakes.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level=loglevel, filemode='a')
    #log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=loglevel)

    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))

logid = 'header.smk: '


# Parse SUBCONFIG
BINS = config["BINS"]
MAXTHREAD = int(config["MAXTHREADS"])
SAMPLES = [os.path.join(x) for x in sampleslong(config)]

if len(SAMPLES) < 1:
    log.error(logid+'No samples found, please check config file')
    sys.exit(logid+'ERROR: No samples found, please check config file')

SETUP = keysets_from_dict(config['SETTINGS'], 'SAMPLES')[0]
SETS = os.sep.join(SETUP)
SETTINGS = subDict(config['SETTINGS'], SETUP)


# Parse SETTINGS
SEQUENCING = SETTINGS.get('SEQUENCING')
REFERENCE = SETTINGS.get('REFERENCE')
REFDIR = str(os.path.dirname(REFERENCE))
INDEX = SETTINGS.get('INDEX')
PREFIX = SETTINGS.get('PREFIX')
ANNO = SETTINGS.get('ANNOTATION')
IP = SETTINGS.get('IP')
rundedup = True if (config.get('RUNDEDUP')) == 'enabled' else False
if rundedup:
    log.debug('DEDUPLICATION ENABLED')

log.info(logid+'Working on SAMPLES: '+str(SAMPLES))

paired = checkpaired([SAMPLES[0]], config)
if paired == 'paired':
    log.info('RUNNING SNAKEMAKE IN PAIRED READ MODE')

stranded = checkstranded([SAMPLES[0]], config)
if stranded != '':
    log.info('RUNNING SNAKEMAKE WITH STRANDEDNESS '+str(stranded))

# MAPPING Variables
if 'MAPPING' in config:
    MAPCONF = subDict(config['MAPPING'], SETUP)
    log.debug(logid+'MAPPINGCONFIG: '+str(SETUP)+'\t'+str(MAPCONF))
    REF = MAPCONF.get('REFERENCE') if MAPCONF.get('REFERENCE') else MAPCONF[XENV].get('REFERENCE')
    MANNO = MAPCONF.get('ANNOTATION') if MAPCONF.get('ANNOTATION') else MAPCONF[XENV].get('ANNOTATION')
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    if MANNO:
        ANNOTATION = MANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO else ANNO.get('GFF')  # by default GTF format will be used
    MAPPERBIN, MAPPERENV = env_bin_from_config3(config, 'MAPPING')
    IDX = MAPCONF.get('INDEX')
    if IDX:
        INDEX = IDX
    if not INDEX:
        INDEX = str.join(os.sep, [REFDIR, 'INDICES', MAPPERENV])+'.idx'
    UIDX = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'][0]))
    INDICES = INDEX.split(',') if INDEX else list(UIDX)
    INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
    PRE = MAPCONF.get('PREFIX')
    if PRE:
        PREFIX = PRE
    if not PREFIX:
        PREFIX = MAPPERENV
    if len(INDICES) > 1:
        if str(os.path.abspath(INDICES[1])) not in UIDX:
            INDEX2 = str(os.path.abspath(INDICES[1]))
        else:
            INDEX2 = str(os.path.abspath(INDICES[1]))+'_idx'
    else:
        INDEX2 = ''


# Peak Calling Variables
if 'PEAKS' in config:
    PEAKCONF = subDict(config['PEAKS'], SETUP)
    REF = PEAKCONF.get('REFERENCE') if PEAKCONF.get('REFERENCE') else PEAKCONF[XENV].get('REFERENCE')
    ANNOPEAK = PEAKCONF.get('ANNOTATION') if PEAKCONF.get('ANNOTATION') else PEAKCONF[XENV].get('ANNOTATION')
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    if ANNOPEAK:
        ANNOTATION = ANNOPEAK
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO else ANNO.get('GFF')  # by default GTF forma
    if not IP:
        IP = check_ip(SAMPLES, config)
    log.info(logid+'Running Peak finding for '+IP+' protocol')


# UCSC/COUNTING Variables
for x in ['UCSC', 'COUNTING']:
    if x in config:
        XBIN, XENV = env_bin_from_config3(config, x)
        XCONF = subDict(config[x], SETUP)
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF.get('REFERENCE') if XCONF.get('REFERENCE') else XCONF[XENV].get('REFERENCE')
        XANNO = XCONF.get('ANNOTATION') if XCONF.get('ANNOTATION') else XCONF[XENV].get('ANNOTATION')
        if XANNO:
            ANNOTATION = XANNO
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO else ANNO.get('GFF')  # by default GTF forma
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if XENV == 'salmon':
            IDX = XCONF.get('INDEX')
            if IDX:
                INDEX = IDX
            if not INDEX:
                INDEX = str.join(os.sep, [REFDIR, 'INDICES', XENV])+'.idx'
                UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe=XENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'][0]))
            INDICES = INDEX.split(',') if INDEX else list(UIDX)
            INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'


# DE/DEU/DAS/DTU Variables
for x in ['DE', 'DEU', 'DAS', 'DTU']:
    if x in config:
        XCONF = subDict(config[x], SETUP)
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF.get('REFERENCE') if XCONF.get('REFERENCE') else XCONF[XENV].get('REFERENCE')
        XANNO = XCONF.get('ANNOTATION') if XCONF.get('ANNOTATION') else XCONF[XENV].get('ANNOTATION')
        if XANNO:
            ANNOTATION = XANNO
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO else ANNO.get('GFF')  # by default GTF forma
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))


# CIRCS Variables
if 'CIRCS' in config:
    CIRCCONF = subDict(config['CIRCS'], SETUP)
    log.debug(logid+'CIRCCONFIG: '+str(SETUP)+'\t'+str(CIRCCONF))
    REF = CIRCCONF.get('REFERENCE') if CIRCCONF.get('REFERENCE') else CIRCCONF[XENV].get('REFERENCE')
    XANNO = CIRCCONF.get('ANNOTATION') if CIRCCONF.get('ANNOTATION') else CIRCCONF[XENV].get('ANNOTATION')
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    CANNO = CIRCCONF.get('ANNOTATION')
    if CANNO:
        ANNOTATION = CANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO else ANNO.get('GFF')  # by default GTF format will be used


combo = ''

####HEADER ENDS HERE####
