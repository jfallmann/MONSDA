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

from MONSDA.Logger import makelogdir, setup_logger
from MONSDA.Params import checkpaired, checkstranded, check_IP, sampleslong, basecall_samples, download_samples, env_bin_from_config, tool_params, samplecond, conditiononly, comparable_as_string, get_cutoff_as_string, get_reps, get_diego_samples, get_diego_groups, get_pairing, set_pairing, get_samples_postprocess, fixRunParameters, resolve_sequencing
from MONSDA.Utils import get_dict_hash, sub_dict, makeoutdir, keysets_from_dict, dict_inst

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

#Make sure we have a TMP directory
makeoutdir("TMP")

# Parse SUBCONFIG
try:
    installpath = os.path.dirname(__file__).replace(
        os.sep.join(["lib", "python3.9", "site-packages", "MONSDA"]), "share"
    )
except:
    installpath = os.path.cwd()

BINS = config.get("BINS")
MAXTHREAD = int(config["MAXTHREADS"])

# set maximum available memory for sorting
def get_sortmem(w, resources):
    # Otherwise, try mem of form 8G / 8GB, again skipping `<TBD>`
    if hasattr(resources, "mem"):
        mem_val = str(resources.mem)
        if mem_val != "<TBD>":
            mem_str = re.sub(r"GB?$", "", mem_val, flags=re.IGNORECASE).strip()
            return max(int(mem_str) -2, 4)  # reserve 2GB for other processes, min 4GB

    # Convert mem_mb if present and not `<TBD>`
    if hasattr(resources, "mem_mb"):
        mem_val = str(resources.mem_mb)
        if mem_val != "<TBD>":
            mem_mb = int(mem_val)
            # round up to full GB
            return max(int((mem_mb + 1023) / 1024) - 2, 4)  # reserve 2GB for other processes, min 4GB
    
    # Fallback if everything is `<TBD>` or missing
    # pick something conservative but nonzero
    return 4  # 4 GB

if not config.get('FETCH', False) and not config.get("BASECALL", False):
    SAMPLES = [os.path.join(x) for x in sampleslong(config)]  
elif not config.get("BASECALL", False):
    SAMPLES = [os.path.join(x) for x in download_samples(config)] 
else:
    SAMPLES = [os.path.join(x) for x in basecall_samples(config)]

if len(SAMPLES) < 1:
    log.error(logid+'No samples found, please check config file')
    sys.exit(logid+'ERROR: No samples found, please check config file')

SETUP = keysets_from_dict(config['SETTINGS'], 'SAMPLES')[0]
SETS = os.sep.join(SETUP)
SETTINGS = sub_dict(config['SETTINGS'], SETUP)


# Parse SETTINGS
SEQUENCING = SETTINGS.get('SEQUENCING')
REFERENCE = SETTINGS.get('REFERENCE')
REFDIR = str(os.path.dirname(REFERENCE))
INDEX = SETTINGS.get('INDEX')
INDEX2 = SETTINGS.get('INDEX2')
PREFIX = SETTINGS.get('PREFIX')
ANNO = SETTINGS.get('ANNOTATION')
DECOY = SETTINGS.get("DECOY")
IP = SETTINGS.get('IP')
rundedup = True if (config.get('RUNDEDUP')) == 'enabled' else False
prededup = True if (config.get('PREDEDUP')) == 'enabled' else False

ANNOTATION = False  # by default no annotation is set for every workflow step

if rundedup:
    if prededup:
        log.info(logid+'(PRE)DEDUPLICATION ENABLED')
    else:
        log.info(logid+'DEDUPLICATION ENABLED')

log.info(logid+'Working on SAMPLES: '+str(SAMPLES))

paired = checkpaired([SAMPLES[0]], config)
if paired == 'paired':
    log.info('RUNNING SNAKEMAKE IN PAIRED READ MODE')
elif paired == 'singlecell':
    log.info('RUNNING SNAKEMAKE IN Singlecell MODE')
stranded = checkstranded([SAMPLES[0]], config)
if stranded != '':
    log.info('RUNNING SNAKEMAKE WITH STRANDEDNESS '+str(stranded))

def raw_fastq_input(wildcards, read=None):
    sample = getattr(wildcards, 'file', None) or getattr(wildcards, 'rawfile')
    rawfile = [x for x in SAMPLES if x.split(os.sep)[-1] in sample][0]
    read = read or getattr(wildcards, 'read', None)
    suffix = f"_{read}" if read else ("_R1" if paired in ['paired', 'singlecell'] else "")
    return f"FASTQ/{rawfile}{suffix}.fastq.gz"

def dedup_fastq_input(wildcards, read=None):
    read = read or getattr(wildcards, 'read', None)
    suffix = f"_{read}" if read else ("_R1" if DEDUP_PAIRED == 'paired' else "")
    return f"DEDUP_FASTQ/{wildcards.combo}/{wildcards.file}{suffix}_dedup.fastq.gz"

def trimmed_fastq_input(wildcards, read=None):
    read = read or getattr(wildcards, 'read', None)
    suffix = f"_{read}" if read else ("_R1" if TRIMMING_PAIRED == 'paired' else "")
    return f"TRIMMED_FASTQ/{wildcards.combo}/{wildcards.file}{suffix}_trimmed.fastq.gz"


# Per-step pairedness/strandedness, default to the global values above and
# overridden below wherever a workflow step defines its own SEQUENCING key.
QC_PAIRED, QC_STRANDED = paired, stranded
TRIMMING_PAIRED, TRIMMING_STRANDED = paired, stranded
DEDUP_PAIRED, DEDUP_STRANDED = paired, stranded
MAPPING_PAIRED, MAPPING_STRANDED = paired, stranded
PEAKS_PAIRED, PEAKS_STRANDED = paired, stranded
TRACKS_PAIRED, TRACKS_STRANDED = paired, stranded
COUNTING_PAIRED, COUNTING_STRANDED = paired, stranded
CIRCS_PAIRED, CIRCS_STRANDED = paired, stranded
FUSIONS_PAIRED, FUSIONS_STRANDED = paired, stranded
DE_PAIRED, DE_STRANDED = paired, stranded
DEU_PAIRED, DEU_STRANDED = paired, stranded
DAS_PAIRED, DAS_STRANDED = paired, stranded
DTU_PAIRED, DTU_STRANDED = paired, stranded

for step in ['QC', 'TRIMMING', 'DEDUP']:
    if step not in config:
        continue
    if step == 'QC':
        stepenv = config.get('POSTQCENV') or config.get('PREQCENV')
    else:
        _, stepenv = env_bin_from_config(config, step)
    if not stepenv:
        continue
    stepenv = stepenv.split('_')[0]
    stepconf = sub_dict(config[step], SETUP)
    step_toolconf = stepconf.get(stepenv, {})
    step_seq = step_toolconf.get('SEQUENCING')
    if not step_seq:
        step_seq = stepconf.get('SEQUENCING') or config[step].get('SEQUENCING')
    globals()[step+'_PAIRED'], globals()[step+'_STRANDED'] = resolve_sequencing(
        step_seq, paired, stranded
    )

# MAPPING Variables
if 'MAPPING' in config:
    MAPPERBIN, MAPPERENV = env_bin_from_config(config, 'MAPPING')
    MAPCONF = sub_dict(config['MAPPING'], SETUP)
    MAPPERENV = MAPPERENV.split('_')[0]
    log.debug(logid+'MAPPINGCONFIG: '+str(SETUP)+'\t'+str(MAPCONF))
    REF = MAPCONF[MAPPERENV].get("REFERENCE", MAPCONF.get("REFERENCE", config['MAPPING'].get("REFERENCE")))
    MANNO = MAPCONF[MAPPERENV].get("ANNOTATION", MAPCONF.get("ANNOTATION", config['MAPPING'].get("ANNOTATION")))
    MDECOY = MAPCONF[MAPPERENV].get("DECOY", MAPCONF.get("DECOY", config['MAPPING'].get("DECOY")))
    MSEQ = MAPCONF[MAPPERENV].get("SEQUENCING", MAPCONF.get("SEQUENCING", config['MAPPING'].get("SEQUENCING")))
    MAPPING_PAIRED, MAPPING_STRANDED = resolve_sequencing(MSEQ, paired, stranded)
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    if MANNO and MANNO != '':
        ANNOTATION = MANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
    if MDECOY and not dict_inst(MDECOY):
            DECOY = os.path.abspath(MDECOY)
    elif dict_inst(DECOY) and DECOY.get(MAPPERENV):
            DECOY = os.path.abspath(DECOY.get(MAPPERENV))
    else:
        DECOY = None
    IDX = MAPCONF.get('INDEX', MAPCONF[MAPPERENV].get('INDEX'))
    if IDX:
        INDEX = IDX
    if not INDEX:
        INDEX = str.join(os.sep, [REFDIR, 'INDICES', MAPPERENV])+'.idx'
    keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
    keydict["REF"] = REFERENCE
    keydict["DECOY"] = DECOY
    keydict["ENV"] = MAPPERENV
    unik = get_dict_hash(keydict)
    UIDX = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=unik)
    INDICES = INDEX.split(',') if INDEX else list(UIDX)
    INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
    IDX2 = MAPCONF.get('INDEX2', MAPCONF[MAPPERENV].get('INDEX2'))
    if IDX2:
        INDEX2 = IDX2
    if not INDEX2 and ('segemehl' in MAPPERENV and 'bisulfite' in MAPPERENV):
        if len(INDICES) > 1:
            UIDX2 = expand("{refd}/INDICES/{mape}/{unikey}.idx2", refd=REFDIR, mape=MAPPERENV, unikey=unik)
            INDEX2 = str(os.path.abspath(INDICES[1])) if str(os.path.abspath(INDICES[1])) not in UIDX2 else str(os.path.abspath(INDICES[1]))+'_idx2'
        else:
            INDEX2 = str.join(os.sep, [REFDIR, 'INDICES', MAPPERENV])+'.idx2'
    else:
        INDEX2 = None
        UIDX2 = None
    MAPOPT = MAPCONF.get(MAPPERENV).get('OPTIONS')
    PRE = MAPCONF.get('PREFIX', MAPCONF.get('EXTENSION', MAPOPT.get('PREFIX', MAPOPT.get('EXTENSION'))))
    if PRE and PRE is not None:
        PREFIX = PRE
    if not PREFIX or PREFIX is None:
        PREFIX = MAPPERENV


# Peak Calling Variables
if 'PEAKS' in config:
    PEAKCONF = sub_dict(config['PEAKS'], SETUP)
    PEAKBIN, PEAKENV = env_bin_from_config(config, 'PEAKS')
    REF = PEAKCONF.get('REFERENCE', PEAKCONF[PEAKENV].get('REFERENCE', config['PEAKS'].get('REFERENCE')))
    ANNOPEAK = PEAKCONF.get('ANNOTATION', PEAKCONF[PEAKENV].get('ANNOTATION', config['PEAKS'].get('ANNOTATION')))
    PSEQ = PEAKCONF.get('SEQUENCING', PEAKCONF[PEAKENV].get('SEQUENCING', config['PEAKS'].get('SEQUENCING')))
    # downstream steps consume the MAPPING output, so default to its pairedness
    PEAKS_PAIRED, PEAKS_STRANDED = resolve_sequencing(PSEQ, MAPPING_PAIRED, MAPPING_STRANDED)
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
        XBIN, XENV = env_bin_from_config(config, x)
        XCONF = sub_dict(config[x], SETUP)
        XENV = XENV.split('_')[0]
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF[XENV].get('REFERENCE', XCONF.get('REFERENCE', config[x].get('REFERENCE')))
        XANNO = XCONF[XENV].get('ANNOTATION', XCONF.get('ANNOTATION', config[x].get('ANNOTATION')))
        XDECOY = XCONF[XENV].get("DECOY", XCONF.get("DECOY", config[x].get("DECOY")))
        XSEQ = XCONF[XENV].get("SEQUENCING", XCONF.get("SEQUENCING", config[x].get("SEQUENCING")))
        # TRACKS/COUNTING consume the MAPPING output, so default to its pairedness
        globals()[x+'_PAIRED'], globals()[x+'_STRANDED'] = resolve_sequencing(XSEQ, MAPPING_PAIRED, MAPPING_STRANDED)
        if XANNO and XANNO != '':
            ANNOTATION = XANNO
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
        if XDECOY and not dict_inst(XDECOY):
            DECOY = XDECOY
        elif dict_inst(DECOY) and DECOY.get(XENV):
            DECOY = DECOY.get(XENV)
        else:
            DECOY = None
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if XENV in ["salmon", "kallisto"]:
            IDX = XCONF.get('INDEX')
            if IDX:
                INDEX = IDX
            if not INDEX:
                INDEX = str.join(os.sep, [REFDIR, 'INDICES', XENV])+'.idx'
                keydict = sub_dict(tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'], ['INDEX'])
                keydict["REF"] = REFERENCE
                keydict["DECOY"] = DECOY
                keydict["ENV"] = XENV
                unik = get_dict_hash(keydict)
                UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe=XENV, unikey=unik)
            INDICES = INDEX.split(',') if INDEX else list(UIDX)
            INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'


# DE/DEU/DAS/DTU Variables
for x in ['DE', 'DEU', 'DAS', 'DTU']:
    if x in config:
        XCONF = sub_dict(config[x], SETUP)
        XBIN, XENV = env_bin_from_config(config, x)
        XENV = XENV.split('_')[0]
        log.debug(logid+'XCONFIG: '+str(SETUP)+'\t'+str(XCONF))
        REF = XCONF[XENV].get("REFERENCE", XCONF.get("REFERENCE", config[x].get("REFERENCE")))
        XANNO = XCONF[XENV].get("ANNOTATION", XCONF.get("ANNOTATION", config[x].get("ANNOTATION")))
        XDECOY = XCONF[XENV].get("DECOY", XCONF.get("DECOY", config[x].get("DECOY")))
        XSEQ = XCONF[XENV].get("SEQUENCING", XCONF.get("SEQUENCING", config[x].get("SEQUENCING")))
        if x == 'DTU':
            # DTU's default salmon path quantifies directly from TRIMMED_FASTQ,
            # independent of the MAPPING step, so fall back to the global layout
            globals()[x+'_PAIRED'], globals()[x+'_STRANDED'] = resolve_sequencing(XSEQ, paired, stranded)
        else:
            globals()[x+'_PAIRED'], globals()[x+'_STRANDED'] = resolve_sequencing(XSEQ, MAPPING_PAIRED, MAPPING_STRANDED)
        if XANNO and XANNO != '':
            ANNOTATION = XANNO            
        else:
            ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
        if XDECOY and not dict_inst(XDECOY):
            DECOY = XDECOY
        elif dict_inst(DECOY) and DECOY.get(XENV):
            DECOY = DECOY.get(XENV)
        else:
            DECOY = None
        if REF:
            REFERENCE = REF
            REFDIR = str(os.path.dirname(REFERENCE))
        if x == 'DTU':
            COUNTTOOL = tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'].get('COUNTTOOL', 'salmon')
            if COUNTTOOL == 'salmon':
                IDX = XCONF.get('INDEX')
                if IDX:
                    INDEX = IDX
                if not INDEX:
                    INDEX = str.join(os.sep, [REFDIR, 'INDICES', XENV])+'.idx'
                    keydict = sub_dict(tool_params(SAMPLES[0], None, config, x, XENV)['OPTIONS'], ['INDEX'])
                    keydict["REF"] = REFERENCE
                    keydict["DECOY"] = DECOY
                    keydict["ENV"] = XENV
                    unik = get_dict_hash(keydict)
                    UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe=COUNTTOOL, unikey=unik)
                INDICES = INDEX.split(',') if INDEX else list(UIDX)
                INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
            else:
                # Long-read quantifiers (e.g. oarfish) quantify directly from
                # the genome-mapped BAM, no transcript index is built.
                INDEX = None
                UIDX = None

# CIRCS Variables
if 'CIRCS' in config:
    CIRCCONF = sub_dict(config['CIRCS'], SETUP)
    XBIN, XENV = env_bin_from_config(config, 'CIRCS')
    log.debug(logid+'CIRCCONFIG: '+str(SETUP)+'\t'+str(CIRCCONF))
    REF = CIRCCONF.get('REFERENCE', CIRCCONF[XENV].get('REFERENCE', config['CIRCS'].get('REFERENCE')))
    XANNO = CIRCCONF.get('ANNOTATION', CIRCCONF[XENV].get('ANNOTATION', config['CIRCS'].get('ANNOTATION')))
    CSEQ = CIRCCONF.get('SEQUENCING', CIRCCONF[XENV].get('SEQUENCING', config['CIRCS'].get('SEQUENCING')))
    CIRCS_PAIRED, CIRCS_STRANDED = resolve_sequencing(CSEQ, MAPPING_PAIRED, MAPPING_STRANDED)
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    CANNO = CIRCCONF.get('ANNOTATION')
    if CANNO and CANNO != '':
        ANNOTATION = CANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used

# FUSIONS Variables
if 'FUSIONS' in config:
    FUSCONF = sub_dict(config['FUSIONS'], SETUP)
    XBIN, XENV = env_bin_from_config(config, 'FUSIONS')
    log.debug(logid+'FUSCONFIG: '+str(SETUP)+'\t'+str(FUSCONF))
    REF = FUSCONF.get('REFERENCE', FUSCONF[XENV].get('REFERENCE', config['FUSIONS'].get('REFERENCE')))
    XANNO = FUSCONF.get('ANNOTATION', FUSCONF[XENV].get('ANNOTATION', config['FUSIONS'].get('ANNOTATION')))
    FSEQ = FUSCONF.get('SEQUENCING', FUSCONF[XENV].get('SEQUENCING', config['FUSIONS'].get('SEQUENCING')))
    FUSIONS_PAIRED, FUSIONS_STRANDED = resolve_sequencing(FSEQ, MAPPING_PAIRED, MAPPING_STRANDED)
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    FANNO = FUSCONF.get('ANNOTATION')
    if FANNO and FANNO != '':
        ANNOTATION = FANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used


if ANNOTATION and isinstance(ANNOTATION, str):
    # Resolve once to an absolute path so every rule's shell command can read
    # it reliably regardless of its own working directory (e.g. containers).
    ANNOTATION = os.path.abspath(ANNOTATION)

combo = ''

if ANNOTATION:
    log.info(logid+f'Using COMBO: {combo}, REF: {REFERENCE}, ANNO: {ANNOTATION}, DECOY: {DECOY}, INDEX: {INDEX}, INDEX2: {INDEX2}, PREFIX: {PREFIX}, IP: {IP}, SETUP: {SETUP}, SEQUENCING: {SEQUENCING}, STRANDED: {stranded}, PAIRED: {paired}, RUNDEDUP: {rundedup}, PREDEDUP: {prededup}, SAMPLES: {SAMPLES}, BINS: {BINS}, MAXTHREAD: {MAXTHREAD}')

####HEADER ENDS HERE####
