import glob, os, sys, inspect, snakemake, json, shutil
import tempfile
import traceback as tb
from collections import defaultdict
from itertools import combinations

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../snakes/lib"),os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"snakes/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from Collection import *
from Logger import *

log = setup_logger(name='snakemake', log_file='LOGS/snakemake.log', level='DEBUG')

logid = 'header.smk: '
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
MAXTHREAD=int(config["MAXTHREADS"])
SOURCE=sources(config)

SAMPLES = [os.path.join(x) for x in sampleslong(config)]
check = [os.path.join('FASTQ',str(x)+'*.fastq.gz') for x in SAMPLES]
SAMPLES = list()
for s in check:
    log.debug(logid+'SEARCHING: '+s)
    f = glob.glob(s)
    log.debug(logid+'SAMPLECHECK: '+str(f))
    if f:
        SAMPLES.extend([str.join(os.sep,x.split(os.sep)[1:]).replace('.fastq.gz','') for x in f])
log.debug(logid+'SAMPLETEST: '+str(SAMPLES))
if len(SAMPLES) < 1:
    log.error(logid+'No samples found, please check config file')
    sys.exit()

log.info(logid+'Working on SAMPLES: '+str(SAMPLES))

paired = checkpaired(SAMPLES, config)
if paired != '':
    log.info('RUNNING SNAKEMAKE IN PAIRED READ MODE')

stranded = checkstranded(SAMPLES, config)
if stranded != '':
    log.info('RUNNING SNAKEMAKE WITH STRANDEDNESS '+str(stranded))

if 'PEAKS' in config:
    CLIP = checkclip(SAMPLES, config)
    peakconf = tool_params(SAMPLES[0],None,config,'PEAKS')['OPTIONS'][0]
    if 'ANNOTATION' in tool_params(SAMPLES[0],None,config,'PEAKS'):
        ANNOPEAK = tool_params(SAMPLES[0],None,config,'PEAKS')['ANNOTATION']
    else:
        ANNOPEAK = None
    try:
        all([x in peakconf for x in ['MINPEAKRATIO', 'PEAKDISTANCE', 'PEAKWIDTH', 'PEAKCUTOFF', 'MINPEAKHEIGHT', 'USRLIMIT']])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error('Not all required options defined in config!\n'+''.join(tbe.format()))

    MINPEAKRATIO = peakconf['MINPEAKRATIO']
    PEAKDISTANCE = peakconf['PEAKDISTANCE']
    PEAKWIDTH = peakconf['PEAKWIDTH']
    PEAKCUTOFF = peakconf['PEAKCUTOFF']
    MINPEAKHEIGHT = peakconf['MINPEAKHEIGHT']
    USRLIMIT = peakconf['USRLIMIT']

    if 'PREPROCESS' in peakconf:
        PREPROCESS = ' '.join("{!s} {!s}".format(key,val) for (key,val) in peakconf['PREPROCESS'].items())
    else:
        PREPROCESS = ''
