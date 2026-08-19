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
from MONSDA.Params import checkpaired, checkstranded, check_IP, sampleslong, basecall_samples, download_samples, env_bin_from_config, tool_params, samplecond, conditiononly, comparable_as_string, get_cutoff_as_string, get_reps, get_diego_samples, get_diego_groups, get_pairing, set_pairing, get_samples_postprocess, fixRunParameters
from MONSDA.Utils import get_dict_hash, sub_dict, makeoutdir, keysets_from_dict, dict_inst
loglevel='INFO'
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
# MAPPING Variables
if 'MAPPING' in config:
    MAPPERBIN, MAPPERENV = env_bin_from_config(config, 'MAPPING')
    MAPCONF = sub_dict(config['MAPPING'], SETUP)
    MAPPERENV = MAPPERENV.split('_')[0]
    log.debug(logid+'MAPPINGCONFIG: '+str(SETUP)+'\t'+str(MAPCONF))
    REF = MAPCONF[MAPPERENV].get("REFERENCE", MAPCONF.get("REFERENCE", config['MAPPING'].get("REFERENCE")))
    MANNO = MAPCONF[MAPPERENV].get("ANNOTATION", MAPCONF.get("ANNOTATION", config['MAPPING'].get("ANNOTATION")))
    MDECOY = MAPCONF[MAPPERENV].get("DECOY", MAPCONF.get("DECOY", config['MAPPING'].get("DECOY")))
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
                UIDX = expand("{refd}/INDICES/{xe}/{unikey}.idx", refd=REFDIR, xe="salmon", unikey=unik)
            INDICES = INDEX.split(',') if INDEX else list(UIDX)
            INDEX = str(os.path.abspath(INDICES[0])) if str(os.path.abspath(INDICES[0])) not in UIDX else str(os.path.abspath(INDICES[0]))+'_idx'
# CIRCS Variables
if 'CIRCS' in config:
    CIRCCONF = sub_dict(config['CIRCS'], SETUP)
    XBIN, XENV = env_bin_from_config(config, 'CIRCS')
    log.debug(logid+'CIRCCONFIG: '+str(SETUP)+'\t'+str(CIRCCONF))
    REF = CIRCCONF.get('REFERENCE', CIRCCONF[XENV].get('REFERENCE', config['CIRCS'].get('REFERENCE')))
    XANNO = CIRCCONF.get('ANNOTATION', CIRCCONF[XENV].get('ANNOTATION', config['CIRCS'].get('ANNOTATION')))
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
    if REF:
        REFERENCE = REF
        REFDIR = str(os.path.dirname(REFERENCE))
    FANNO = FUSCONF.get('ANNOTATION')
    if FANNO and FANNO != '':
        ANNOTATION = FANNO
    else:
        ANNOTATION = ANNO.get('GTF') if 'GTF' in ANNO and ANNO.get('GTF') != '' else ANNO.get('GFF')  # by default GTF format will be used
combo = ''
if ANNOTATION:
    log.info(logid+f'Using COMBO: {combo}, REF: {REFERENCE}, ANNO: {ANNOTATION}, DECOY: {DECOY}, INDEX: {INDEX}, INDEX2: {INDEX2}, PREFIX: {PREFIX}, IP: {IP}, SETUP: {SETUP}, SEQUENCING: {SEQUENCING}, STRANDED: {stranded}, PAIRED: {paired}, RUNDEDUP: {rundedup}, PREDEDUP: {prededup}, SAMPLES: {SAMPLES}, BINS: {BINS}, MAXTHREAD: {MAXTHREAD}')
####HEADER ENDS HERE####

combo = 'fastqc-cutadapt-picard-rustar'

wildcard_constraints:
    combo = combo,
    rawfile = '|'.join(list(SAMPLES)),
    file = '|'.join(list(samplecond(SAMPLES, config))),
    read = "R1|R2"


rule themall:
	input:	expand("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", combo=combo, file=samplecond(SAMPLES, config), type=["sorted", "sorted_unique"]),
		expand("QC/Multi/{combo}/{condition}/multiqc_report.html", condition=str.join(os.sep, conditiononly(SAMPLES[0], config)), combo=combo)

QCBIN, QCENV = env_bin_from_config(config, 'PREQC')
if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_{read}_raw.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_{read}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/{read}_fastqc_dedup.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/{read}_fastqc_trimmed.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_raw.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_dedup.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_trimmed.log"
        conda: "<REPO>/envs/"+QCENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"
rule qc_mapped:
    input:   r1 = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output:  o1 = report("QC/{combo}/{file}_mapped_sorted_fastqc.zip", category="QC")
    log:     "LOGS/{combo}/{file}/QC/fastqc/fastqc_mapped.log"
    conda: "<REPO>/envs/"+QCENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
rule qc_uniquemapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_uniquemapped.log"
    conda: "<REPO>/envs/"+QCENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
rule qc_dedupmapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_dedupmapped.log"
    conda: "<REPO>/envs/"+QCENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
rule qc_uniquededup:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_uniquededup.log"
    conda: "<REPO>/envs/"+QCENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+QCENV+""+"-VERSION"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"


TRIMBIN, TRIMENV = env_bin_from_config(config,'TRIMMING')
if paired == 'paired':
    rule cutadapt_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz"
        output: o1 = "TRIMMED_FASTQ/{combo}/{file}_R1_val_1.fq.gz" if not prededup else "TRIMMED_FASTQ/{combo}/{file}_R1_dedup_val_1.fq.gz",
                o2 = "TRIMMED_FASTQ/{combo}/{file}_R2_val_2.fq.gz" if not prededup else "TRIMMED_FASTQ/{combo}/{file}_R2_dedup_val_2.fq.gz"
        log:    "LOGS/{combo}/{file}/TRIMMING/cutadapt/trim.log"
        conda: "<REPO>/envs/"+TRIMENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+TRIMENV+""+"-VERSION"
        threads: min(int(MAXTHREAD/8),4) if min(int(MAXTHREAD/8),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
        params: odir=lambda wildcards, output:os.path.dirname(output.o1),
                tpara = lambda wildcards: tool_params(wildcards.file, None, config, "TRIMMING", TRIMENV)['OPTIONS'].get('TRIM', ""),
                trim=TRIMBIN
        shell:  "{params.trim} {params.tpara} --cores {threads} -o {output.o1} -p {output.o2} {input.r1} {input.r2} &> {log}"
    rule cutadapt_rename:
        input:  o1 = rules.cutadapt_trim.output.o1,
                o2 = rules.cutadapt_trim.output.o2
        output: r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz"
        conda: "<REPO>/envs/"+TRIMENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+TRIMENV+""+"-VERSION"
        threads: 1
        shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"
else:
    rule cutadapt_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: o1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fq.gz" if not prededup else "TRIMMED_FASTQ/{combo}/{file}_dedup_trimmed.fq.gz"
        log:    "LOGS/{combo}/{file}/TRIMMING/cutadapt/trim.log"
        conda: "<REPO>/envs/"+TRIMENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+TRIMENV+""+"-VERSION"
        threads: min(int(MAXTHREAD/8),4) if min(int(MAXTHREAD/2),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
        params: odir=lambda wildcards, output: os.path.dirname(output.o1),
                tpara = lambda wildcards:tool_params(wildcards.file, None, config, "TRIMMING", TRIMENV)['OPTIONS'].get('TRIM', ""),
                trim=TRIMBIN,
        shell:  "{params.trim} {params.tpara} --cores {threads} -o {output.o1} {input.r1} > {log}"
    rule cutadapt_rename:
        input:  o1 = rules.cutadapt_trim.output.o1
        output: r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"
        conda: "<REPO>/envs/"+TRIMENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+TRIMENV+""+"-VERSION"
        threads: 1
        shell:  "mv {input.o1} {output.r1}"


DEDUPBIN, DEDUPENV = env_bin_from_config(config, 'DEDUP')
rule dedupbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam"
    output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
            bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
            met = report("MAPPED/{combo}/{file}_mapped_{type}_dedup_metrics.txt", category="DEDUP"),
            td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
    log:    "LOGS/{combo}/{file}/DEDUP/picard/dedupbam_{type}.log"
    conda: "<REPO>/envs/"+DEDUPENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+DEDUPENV+""+"-VERSION"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('JAVA', ""),
            dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && {params.dedup} {params.jpara} MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR {output.td} --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.met} {params.dpara} 2> {log} && samtools index {output.bam} 2>> {log} && cp -f {output.met} $(dirname {log})/"

MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)
rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=unik)),
            idxfile = expand("{refd}/INDICES/{mape}_{unikey}/{pref}.idx", refd=REFDIR, mape=MAPPERENV, unikey=unik, pref=PREFIX),
            tmpfa = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV)),
            tmpref = temp(expand("TMP/{mape}/ref.anno", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/MAPPING/{mape}/idx.log", sets=SETS, mape=MAPPERENV)
    conda: "<REPO>/envs/"+MAPPERENV+".yaml"
    container: "oras://ghcr.io/jfallmann/monsda:"+MAPPERENV+""+"-VERSION"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda w,output: fixRunParameters(config, MAPPERENV, SAMPLES[0], None, 'MAPPING', 'INDEX', "--sjdbGTFfile", f"--sjdbGTFfile {output.tmpref}"),
            anno = ANNOTATION,                        
            pref = PREFIX,
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmpfa} && zcat {params.anno} > {output.tmpref} && mkdir -p {output.uidx} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {output.uidx}/{params.pref} --outTmpDir TMP/rustar_generate_index --genomeDir {output.uidx} --genomeFastaFiles {output.tmpfa}  &> {log} && touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx} && cat {output.uidx}/*Log.out >> {log};fi && rm -rf TMP/rustar_generate_index"
if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                log = "LOGS/{combo}/{file}/MAPPING/star/Log.out",
                log_final = "LOGS/{combo}/{file}/MAPPING/star/Log.final.out",
                tmp = temp("TMP/RUSTAROUT/{combo}/{file}")
        log:    "LOGS/{combo}/{file}/MAPPING/star/mapping.log"
        conda: "<REPO>/envs/"+MAPPERENV+".yaml"
        container: "oras://ghcr.io/jfallmann/monsda:"+MAPPERENV+""+"-VERSION"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN,
                anno = ANNOTATION,
                pref = PREFIX,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && gzip -c {output.tmp}.Aligned.out.sam > {output.mapped} && rm -f {output.tmp}.Aligned.out.sam 2>> {log} && touch {output.tmp}.Unmapped.out.mate1 {output.tmp}.Unmapped.out.mate2 && cat {output.tmp}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\" |gzip > {output.unmapped_r1} && cat {output.tmp}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > {output.unmapped_r2} && rm -f {output.tmp}.Unmapped.out* && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp} "
else:
    if paired != 'singlecell':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    uidx = rules.generate_index.output.uidx[0],
                    dummy = rules.generate_index.output.idxfile,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz",
                    log = "LOGS/{combo}/{file}/MAPPING/star/Log.out",
                    log_final = "LOGS/{combo}/{file}/MAPPING/star/Log.final.out",
                    tmp = temp("TMP/RUSTAROUT/{combo}/{file}")                    
            log:    "LOGS/{combo}/{file}/MAPPING/star/mapping.log"
            conda: "<REPO>/envs/"+MAPPERENV+".yaml"
            container: "oras://ghcr.io/jfallmann/monsda:"+MAPPERENV+""+"-VERSION"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),                    
                    mapp = MAPPERBIN,
                    anno = ANNOTATION,
                    pref = PREFIX,
                    tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
            shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && gzip -c {output.tmp}.Aligned.out.sam > {output.mapped} && rm -f {output.tmp}.Aligned.out.sam 2>> {log} touch {output.tmp}.Unmapped.out.mate && gzip {output.tmp}.Unmapped.out.mate* && mv {output.tmp}.Unmapped.out.mate*.gz {output.unmapped} 2>> {log} && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
    else:
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    uidx = rules.generate_index.output.uidx[0],
                    dummy = rules.generate_index.output.idxfile,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                    unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                    unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                    log = "LOGS/{combo}/{file}/MAPPING/star/Log.out",
                    log_final = "LOGS/{combo}/{file}/MAPPING/star/Log.final.out",
                    tmp = temp("TMP/RUSTAROUT/{combo}/{file}")
            log:    "LOGS/{combo}/{file}/MAPPING/star/mapping.log"
            conda: "<REPO>/envs/"+MAPPERENV+".yaml"
            container: "oras://ghcr.io/jfallmann/monsda:"+MAPPERENV+""+"-VERSION"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    stranded = lambda x: '--soloStrand Forward' if stranded == 'fr' else '--soloStrand Reverse' if stranded == 'rf' else '--soloStrand Unstranded',
                    mapp = MAPPERBIN,
                    anno = ANNOTATION,
                    pref = PREFIX,
                    r1 = lambda wildcards, input: input.r1 if not '--soloBarcodeMate 1' in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', "") else input.r2.replace('TRIMMED_','').replace('_trimmed.fastq', '.fastq'),
                    r2 = lambda wildcards, input: input.r2 if not '--soloBarcodeMate 1' in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', "") else input.r1.replace('TRIMMED_','').replace('_trimmed.fastq', '.fastq').replace({combo}, ''),
                    tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
            shell: "{params.mapp} {params.mpara} {params.stranded} --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {params.r1} {params.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && samtools view -h {output.tmp}.Aligned.sortedByCoord.out.bam | gzip > {output.mapped} && rm -f {output.tmp}.Aligned.sortedByCoord.out.bam; touch {output.tmp}.Unmapped.out.mate1 {output.tmp}.Unmapped.out.mate2 && cat {output.tmp}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\" |gzip > {output.unmapped_r1} && cat {output.tmp}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > {output.unmapped_r2} && rm -f {output.tmp}.Unmapped.out* && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp} "


rule sortsam:
    input:  mapps = rules.mapping.output.mapped
    output: sortedsam = report("MAPPED/{combo}/{file}_mapped_sorted.sam.gz", category="SORTING"),
            tmphead = temp("MAPPED/{combo}/{file}_mapped_header.gz"),
            tmpfile = temp("TMP/{combo}/{file}")
    log:    "LOGS/{combo}/{file}/MAPPING/samtools/sortsam.log"
    conda: "<REPO>/envs/samtools.yaml"
    container: "oras://ghcr.io/jfallmann/monsda:samtools"+"-VERSION"
    threads: MAXTHREAD
    priority: 100
    params: linkto = lambda wildcards, output: os.path.basename(output.sortedsam),
            sortmem = get_sortmem
    shell: "set +o pipefail; export LC_ALL=C; samtools view -H {input.mapps}|grep '^@HD' |pigz -p 1 -f > {output.tmphead} 2> {log}; samtools view -H {input.mapps}|grep '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p 1 -f >> {output.tmphead} 2>> {log}; samtools view -H {input.mapps}|grep '^@RG'|pigz -p 1 -f >> {output.tmphead} 2>> {log}; samtools view -H {input.mapps}|grep -P '^@PG'|pigz -p {threads} -f >> {output.tmphead} 2>> {log}; samtools view -h {input.mapps} | grep -v \"^@\"|sort --parallel={threads} -S {params.sortmem}G -T TMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output.tmpfile} 2>> {log}; cat {output.tmphead} {output.tmpfile} > {output.sortedsam} 2>> {log}"
rule sam2bam:
    input:  sortedsam = rules.sortsam.output.sortedsam
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted.bam", category="2BAM"),
            bamindex = "MAPPED/{combo}/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{combo}/{file}/MAPPING/samtools/sam2bam.log"
    conda: "<REPO>/envs/samtools.yaml"
    container: "oras://ghcr.io/jfallmann/monsda:samtools"+"-VERSION"
    threads: MAXTHREAD
    params: bins = BINS
    shell: "zcat {input.sortedsam} | samtools view -bS - > {output.bam} && samtools index {output.bam} 2> {log}"
rule uniqsam:
    input:  sortedsam = rules.sortsam.output.sortedsam,
            bam = rules.sam2bam.output
    output: uniqsam = report("MAPPED/{combo}/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{combo}/{file}/MAPPING/samtools/uniqsam.log"
    conda: "<REPO>/envs/samtools.yaml"
    container: "oras://ghcr.io/jfallmann/monsda:samtools"+"-VERSION"
    threads: MAXTHREAD
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input.sortedsam} {output.uniqsam} {threads} 2> {log}"
rule sam2bamuniq:
    input: uniqsam = rules.uniqsam.output,
           bam = rules.sam2bam.output
    output:  uniqbam = report("MAPPED/{combo}/{file}_mapped_sorted_unique.bam", category="2BAM"),
             uniqbamindex = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{combo}/{file}/MAPPING/samtools/sam2bamuniq.log"
    conda: "<REPO>/envs/samtools.yaml"
    container: "oras://ghcr.io/jfallmann/monsda:samtools"+"-VERSION"
    threads: MAXTHREAD
    priority: 50
    params: bins = BINS
    shell: "zcat {input.uniqsam} | samtools view -bS - > {output.uniqbam} && samtools index {output.uniqbam} 2> {log}"


if rundedup:
    if paired == 'paired':
        if prededup:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                        expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedupmapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquededup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"]),
                        versions = "LOGS/versions.txt"
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
                conda: "<REPO>/envs/qc.yaml"
                container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
        else:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),                    
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedupmapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquededup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"]),
                        versions = "LOGS/versions.txt"
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
                conda: "<REPO>/envs/qc.yaml"
                container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
    else:
        if prededup:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedupmapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquededup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo,type=["sorted", "sorted_unique"]),
                        versions = "LOGS/versions.txt"
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
                conda: "<REPO>/envs/qc.yaml"
                container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"                    
        else:
            rule multiqc:
                input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                        expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),                 
                        expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_dedupmapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.qc_uniquededup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                        expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES, config), combo=combo, type=["sorted", "sorted_unique"]),
                        versions = "LOGS/versions.txt"
                output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                        tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                        lst = "QC/Multi/{combo}/{condition}/qclist.txt"
                log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
                conda: "<REPO>/envs/qc.yaml"
                container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
                threads: 1
                params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
                shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
else:
    if paired == 'paired':
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                    versions = "LOGS/versions.txt"
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
            conda: "<REPO>/envs/qc.yaml"
            container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"
    else:
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES, config), combo=combo),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES, config), combo=combo),
                    versions = "LOGS/versions.txt"
            output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                    lst = "QC/Multi/{combo}/{condition}/qclist.txt"
            log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc.log"
            conda: "<REPO>/envs/qc.yaml"
            container: "oras://ghcr.io/jfallmann/monsda:qc"+"-VERSION"
            threads: 1
            params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
            shell:  "OUT=$(dirname {output.html}); SCAN=LOGS/{wildcards.combo}/{wildcards.condition}; COLLECT=$SCAN/MULTIQC/collect; mkdir -p $COLLECT $OUT; for i in {input};do case \"${{i}}\" in *versions.txt) continue;; esac; echo $(dirname \"${{i}}\") >> {output.tmp}; ln -sfnr \"${{i}}\" $COLLECT/ 2>> {log} || cp -f \"${{i}}\" $COLLECT/ 2>> {log};done; cat {output.tmp} |sort -u > {output.lst}; MODS=$(grep -v '^#' {input.versions}|cut -f3|tr ',' '\\n'|grep -vx '-'|sort -u|sed 's/^/-m /'|tr '\\n' ' ');export LC_ALL=C.UTF-8; multiqc -f {params.qpara} $MODS -k json -z -s -o $OUT $SCAN 2>> {log}; cp -f {input.versions} $OUT/"


## Queue-agnostic fallback resources for all rules.
## These defaults apply when a rule has no explicit resource set.
## They keep local/non-Slurm runs from lacking baseline resource values.
_fallback_rule_resources = {
	"mem_mb": 20000,
	"runtime": 480,
	"nodes": 1,
}
_user_rule_defaults = config.get("DEFAULT_RULE_RESOURCES", {})
if isinstance(_user_rule_defaults, dict):
	_fallback_rule_resources.update(_user_rule_defaults)
def _has_resource(rule_obj, resource_name):
	try:
		return resource_name in rule_obj.resources.keys()
	except Exception:
		return hasattr(rule_obj.resources, resource_name)
def _set_resource(rule_obj, resource_name, value):
	try:
		rule_obj.resources[resource_name] = value
	except Exception:
		setattr(rule_obj.resources, resource_name, value)
for _rule in workflow.rules:
	for _res_name, _res_value in _fallback_rule_resources.items():
		if not _has_resource(_rule, _res_name):
			_set_resource(_rule, _res_name, _res_value)
onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str(log))


