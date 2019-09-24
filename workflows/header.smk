import glob, os, sys, inspect, snakemake, json, shutil
from collections import defaultdict
from snakemake.utils import validate, min_version
min_version("5.5.2")

###snakemake -n -j 20 --use-conda -s Workflow/workflows/mapping_paired.smk
###--configfile Workflow/config_compare.json --directory ${PWD}
###--printshellcmds 2> run.log
cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"snakes/lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Collection import *
from Logger import *

log = setup_logger(name='snakemake', log_file='LOGS/snakemake.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='DEBUG')
#log.debug(cmd_subfolder)
#log.debug(sys.path)

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

paired = ''
if checkpaired(SAMPLES, config):
    paired = 'paired'
