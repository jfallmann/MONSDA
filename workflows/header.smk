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

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../MONSDA/MONSDA"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"MONSDA/MONSDA"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../MONSDA"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"MONSDA")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from MONSDA.Logger import *
from MONSDA.Params import *


loglevel="INFO"

try:
    scriptname = os.path.basename(inspect.stack()[-1].filename).replace('.py','')
    if any(x in scriptname for x in ['MONSDA','Configurator']):
        log = logging.getLogger(scriptname)
    else:
        log = logging.getLogger('snakemake')
        for handler in log.handlers[:]:
            handler.close()
            log.removeHandler(handler)
    handler = logging.FileHandler('LOGS/MONSDA.log', mode='a')
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M'))
    log.addHandler(handler)
    lvl = loglevel
    log.setLevel(lvl)

except Exception as err:
    log = setup_logger(name='MONSDA.header', log_file='LOGS/MONSDA.log', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level=loglevel, filemode='a')
    #log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=loglevel)

    exc_type, exc_value, exc_tb = sys.exc_info()
    tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
    )
    log.error(''.join(tbe.format()))

logid = 'header.smk: '


# Parse SUBCONFIG
try:
    installpath = os.path.dirname(__file__).replace(
        os.sep.join(["lib", "python3.9", "site-packages", "MONSDA"]), "share"
    )
except:
    installpath = os.path.cwd()

BINS = config.get("BINS")
MAXTHREAD = int(config["MAXTHREADS"])
SAMPLES = [os.path.join(x) for x in sampleslong(config)] if not config.get('FETCH', False) else ([os.path.join(x) for x in download_samples(config)] if not config.get("BASECALL", False)
        else [os.path.join(x) for x in basecall_samples(config)])

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
prededup = True if (config.get('PREDEDUP')) == 'enabled' else False

if rundedup:
    if prededup:
        log.info(logid+'(PRE)DEDUPLICATION ENABLED')
    else:
        log.info(logid+'DEDUPLICATION ENABLED')

log.info(logid+'Working on SAMPLES: '+str(SAMPLES))

paired = checkpaired([SAMPLES[0]], config)
if paired == 'paired':
    log.info('RUNNING SNAKEMAKE IN PAIRED READ MODE')

stranded = checkstranded([SAMPLES[0]], config)
if stranded != '':
    log.info('RUNNING SNAKEMAKE WITH STRANDEDNESS '+str(stranded))

# MAPPING Variables
if 'MAPPING' in config:
    MAPPERBIN, MAPPERENV = env_bin_from_config3(config, 'MAPPING')
    MAPCONF = subDict(config['MAPPING'], SETUP)
    log.debug(logid+'MAPPINGCONFIG: '+str(SETUP)+'\t'+str(MAPCONF))
    REF = MAPCONF.get('REFERENCE', MAPCONF[MAPPERENV].get('REFERENCE'))
    MANNO = MAPCONF.get('ANNOTATION', MAPCONF[MAPPERENV].get('ANNOTATION'))
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    if MANNO and MANNO != '':
        ANNOTATION = MANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
    IDX = MAPCONF.get('INDEX', MAPCONF[MAPPERENV].get('INDEX'))
    if IDX:
        INDEX = IDX
    if not INDEX:
        INDEX = str.join(os.sep, [REFDIR, 'INDICES', MAPPERENV])+'.idx'
    UIDX = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])))
    INDICES = INDEX.split(',') if INDEX else list(UIDX)
    INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
    MAPOPT = MAPCONF.get(MAPPERENV).get('OPTIONS')
    PRE = MAPCONF.get('PREFIX', MAPCONF.get('EXTENSION', MAPOPT.get('PREFIX', MAPOPT.get('EXTENSION'))))
    if PRE and PRE is not None:
        PREFIX = PRE
    if not PREFIX or PREFIX is None:
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
    PEAKBIN, PEAKENV = env_bin_from_config3(config, 'PEAKS')
    REF = PEAKCONF.get('REFERENCE', PEAKCONF[PEAKENV].get('REFERENCE'))
    ANNOPEAK = PEAKCONF.get('ANNOTATION', PEAKCONF[PEAKENV].get('ANNOTATION'))
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    if ANNOPEAK and ANNOPEAK != '':
        ANNOTATION = ANNOPEAK
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
    if not IP:
        IP = check_IP(SAMPLES, config)
    log.info(logid+'Running Peak finding for '+IP+' protocol')


# TRACKS/COUNTING Variables
for x in ['TRACKS', 'COUNTING']:
    if x in config:
        XBIN, XENV = env_bin_from_config3(config, x)
        XCONF = subDict(config[x], SETUP)
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF.get('REFERENCE', XCONF[XENV].get('REFERENCE'))
        XANNO = XCONF.get('ANNOTATION', XCONF[XENV].get('ANNOTATION'))
        if XANNO and XANNO != '':
            ANNOTATION = XANNO
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if XENV == 'salmon':
            IDX = XCONF.get('INDEX')
            if IDX:
                INDEX = IDX
            if not INDEX:
                INDEX = str.join(os.sep, [REFDIR, 'INDICES', XENV])+'.idx'
                UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe=XENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'], ['INDEX'])))
            INDICES = INDEX.split(',') if INDEX else list(UIDX)
            INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'


# DE/DEU/DAS/DTU Variables
for x in ['DE', 'DEU', 'DAS', 'DTU']:
    if x in config:
        XCONF = subDict(config[x], SETUP)
        XBIN, XENV = env_bin_from_config3(config, x)
        XENV = XENV.split('_')[0]
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF.get('REFERENCE', XCONF[XENV].get('REFERENCE'))
        XANNO = XCONF.get('ANNOTATION', XCONF[XENV].get('ANNOTATION'))
        if XANNO and XANNO != '':
            ANNOTATION = XANNO
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if x == 'DTU':
            IDX = XCONF.get('INDEX')
            if IDX:
                INDEX = IDX
            if not INDEX:
                INDEX = str.join(os.sep, [REFDIR, 'INDICES', XENV])+'.idx'
                UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe="salmon", unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'], ['INDEX'])))
            INDICES = INDEX.split(',') if INDEX else list(UIDX)
            INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'

# CIRCS Variables
if 'CIRCS' in config:
    CIRCCONF = subDict(config['CIRCS'], SETUP)
    XBIN, XENV = env_bin_from_config3(config, 'CIRCS')
    log.debug(logid+'CIRCCONFIG: '+str(SETUP)+'\t'+str(CIRCCONF))
    REF = CIRCCONF.get('REFERENCE', CIRCCONF[XENV].get('REFERENCE'))
    XANNO = CIRCCONF.get('ANNOTATION', CIRCCONF[XENV].get('ANNOTATION'))
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    CANNO = CIRCCONF.get('ANNOTATION')
    if CANNO and CANNO != '':
        ANNOTATION = CANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used


combo = ''

####HEADER ENDS HERE####
