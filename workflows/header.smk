import glob, os, sys, inspect, snakemake
from snakemake.utils import validate, min_version
min_version("5.3.0")

###snakemake -n -j 20 --use-conda -s Workflow/workflows/mapping_paired.smk
###--configfile Workflow/config_compare.json --directory ${PWD}
###--printshellcmds 2> run.log

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Collection import *

QC=config["QC"]
ADAPTERS=config["ADAPTERS"]
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
MAPPERENV=config["MAPPERENV"]
MAPPERBIN=config["MAPPERBIN"]
TRIMENV=config["TRIMENV"]
TRIMBIN=config["TRIMBIN"]
MAXTHREAD=int(config["MAXTHREADS"])
SOURCE=sources(config)
SAMPLES=list(set(samples(config)))
if os.path.exists(SAMPLES[0]) is False:
    SAMPLES=list(set(sampleslong(config)))
