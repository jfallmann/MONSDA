import glob, os, sys, inspect, snakemake, json, shutil
import traceback as tb
from collections import defaultdict
from itertools import combinations

cmd_subfolder = [os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../snakes/lib"),os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"snakes/lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib"), os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"lib")]
for x in cmd_subfolder:
    if x not in sys.path:
        sys.path.insert(0, x)

from Collection import *
from Logger import *

log = setup_logger(name='snakemake', log_file='LOGS/snakemake.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')

REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
MAXTHREAD=int(config["MAXTHREADS"])
SOURCE=sources(config)
SAMPLES=list(set(samples(config)))
if os.path.exists(SAMPLES[0]) is False:
    SAMPLES=list(set(sampleslong(config)))

log.debug('SAMPLES: '+str(SAMPLES))

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
