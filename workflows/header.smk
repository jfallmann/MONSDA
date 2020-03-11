import glob, os, sys, inspect, snakemake, json, shutil, logging
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

try:
    log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))
    if (log.hasHandlers()):
        logr.handlers.clear()
    log.addHandler(logging.FileHandler('LOGS/RunSnakemake.log'))  # streamlog
    log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
except:
    log = setup_logger(name='Snakemake', log_file='stdout', logformat='%(asctime)s %(levelname)-8s %(name)-12s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')#, filemode='a')
    #log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog

logid = 'header.smk: '
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
MAXTHREAD=int(config["MAXTHREADS"])
SOURCE=sources(config)

SAMPLES = [os.path.join(x) for x in sampleslong(config)]
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
        log.error('Not all required peak finding options defined in config!\n'+''.join(tbe.format()))

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
