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

combo = 'fastqc-trimgalore-picard-star/drimseq_DTU'

scombo = 'fastqc-trimgalore-picard-star'

wildcard_constraints:
    combo = combo,
    scombo = scombo,
    read = "R1|R2",
    type = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"

logid = 'drimseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config(config,'DTU')
log.debug(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config(SAMPLES, config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers
comparison = comparable_as_string(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
usededup = config.get('RUNDEDUP', False)
TERMINUSENV = 'terminus'
termpara = tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('TERMINUS', None)
runterminus = termpara is not None
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = DTUENV
unik = get_dict_hash(keydict)
rule themall:
    input:  session = expand("DTU/{combo}/DTU_DRIMSEQ_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            res_t   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_g   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_p   = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_proportions.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            res_gwp = expand("DTU/{combo}/Tables/DTU_DRIMSEQ_{scombo}_{comparison}_table_genewise-precision.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            fig_F   = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_FeatPerGene.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_P   = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_Precision.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_PV  = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_figure_PValues.png", combo=combo, comparison=compstr, scombo=scombo),
            fig_files = expand("DTU/{combo}/Figures/DTU_DRIMSEQ_{scombo}_{comparison}_list_sig{feature}Figures.tsv", combo=combo, comparison=compstr, scombo=scombo, feature=["Gene","Transcript"]),
            sig_t    = expand("DTU/{combo}/Tables/Sig_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_dt   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_ut   = expand("DTU/{combo}/Tables/SigUP_DTU_DRIMSEQ_{scombo}_{comparison}_table_transcripts.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_g    = expand("DTU/{combo}/Tables/Sig_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_dg   = expand("DTU/{combo}/Tables/SigDOWN_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_ug   = expand("DTU/{combo}/Tables/SigUP_DTU_DRIMSEQ_{scombo}_{comparison}_table_genes.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)
rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=unik))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  "<REPO>/envs/"+COUNTENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('INDEX', ""),
            decoy = f"-d {os.path.abspath(DECOY)}" if DECOY else '',
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} {params.decoy} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"
if paired == 'paired':
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{scombo}/{file}_R2_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1)),
                filetolink2 = lambda w, input: "{r}".format(r=os.path.abspath(input.r2))
        shell:  "ln -s {params.filetolink} {output.r1} && ln -s {params.filetolink2} {output.r2}"
else:
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not usededup else "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{scombo}/{file}_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w, input: "{r}".format(r=os.path.abspath(input.r1))
        shell:  "ln -s {params.filetolink} {output.r1}"
if paired == 'paired':
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo),
                r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "<REPO>/envs/"+COUNTENV+".yaml"
        container: "oras://docker.io/jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else '-l IU',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"
else:
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "<REPO>/envs/"+COUNTENV+".yaml"
        container: "oras://docker.io/jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -r {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"
if runterminus:
    rule terminus_group:
        input:  ctsdir = rules.mapping.output.ctsdir
        output: grp = "DTU/{combo}/terminus/{file}/groups.txt"
        log:    "LOGS/{combo}/{file}/terminus_group.log"
        conda:  "<REPO>/envs/"+TERMINUSENV+".yaml"
        container: "oras://docker.io/jfallmann/monsda:"+TERMINUSENV+""
        threads: MAXTHREAD
        params: tpara = termpara,
                outdir = lambda wildcards, output: os.path.dirname(os.path.dirname(str(output.grp)))
        shell:  "terminus group {params.tpara} -d {input.ctsdir} -o {params.outdir} &> {log}"
    rule terminus_collapse:
        input:  grps = expand("DTU/{combo}/terminus/{file}/groups.txt", combo=combo, file=samplecond(SAMPLES, config)),
                dirs = expand(rules.mapping.output.ctsdir, combo=combo, file=samplecond(SAMPLES, config))
        output: cnts = expand("DTU/{combo}/terminus/{file}/quant.sf.gz", combo=combo, file=samplecond(SAMPLES, config))
        log:    expand("LOGS/DTU/{combo}/terminus_collapse.log", combo=combo)
        conda:  "<REPO>/envs/"+TERMINUSENV+".yaml"
        container: "oras://docker.io/jfallmann/monsda:"+TERMINUSENV+""
        threads: MAXTHREAD
        params: outdir = lambda wildcards, input: os.path.dirname(os.path.dirname(str(input.grps[0])))
        shell:  "terminus collapse -d {input.dirs} -o {params.outdir} &> {log}; for d in {input.dirs}; do b=$(basename $d); gzip -c {params.outdir}/$b/quant.sf > {params.outdir}/$b/quant.sf.gz; done 2>> {log}"
    ctsource = expand("DTU/{combo}/terminus/{file}", combo=combo, file=samplecond(SAMPLES, config))
else:
    ctsource = expand(rules.mapping.output.ctsdir, combo=combo, file=samplecond(SAMPLES, config))
rule create_annotation_table:
    input:  dir  = ctsource,
    output: anno = expand("DTU/{combo}/Tables/{scombo}_ANNOTATION.gz", combo=combo, scombo=scombo)
    log:    expand("LOGS/DTU/{combo}/create_DTU_table.log", combo=combo)
    conda:  "<REPO>/envs/"+COUNTENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+COUNTENV+""
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.dir, config, 'DTU'),
            bins = BINS
    shell:  "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"
rule run_DTU:
    input:  anno = expand(rules.create_annotation_table.output.anno, combo=combo, scombo=scombo)
    output: session = rules.themall.input.session,
            res_t   = rules.themall.input.res_t,
            res_g   = rules.themall.input.res_g,
            res_p   = rules.themall.input.res_p,
            res_gwp = rules.themall.input.res_gwp,
            fig_F   = rules.themall.input.fig_F,
            fig_P   = rules.themall.input.fig_P,
            fig_PV  = rules.themall.input.fig_PV,
            fig_files = rules.themall.input.fig_files
    log:    expand("LOGS/DTU/{combo}/run_DTU.log", combo=combo)
    conda:  "<REPO>/envs/"+DTUENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+DTUENV+""
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            outdir = 'DTU/'+combo,
            ref = os.path.abspath(ANNOTATION),
            dtuopt = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('DTU', "")
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.compare} {params.pcombo} {threads} \'{params.dtuopt}\' 2> {log}"
rule filter_significant_drimseq:
    input:  res_g = rules.themall.input.res_g,
            res_t = rules.themall.input.res_t
    output: sig_g   = rules.themall.input.sig_g,
            sig_dg  = rules.themall.input.sig_dg,
            sig_ug  = rules.themall.input.sig_ug,
            sig_t   = rules.themall.input.sig_t,
            sig_dt  = rules.themall.input.sig_dt,
            sig_ut  = rules.themall.input.sig_ut
    log:    "LOGS/DTU/filter_drimseqDTU.log"
    conda:  "<REPO>/envs/"+DTUENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+DTUENV+""
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DTU', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DTU', 'lfc')
    shell:  "set +o pipefail; arr=({input.res_g}); arrt=({input.res_t}); orr=({output.sig_g}); orrd=({output.sig_dg}); orru=({output.sig_ug}); orrt=({output.sig_t}); orrtd=({output.sig_dt}); orrtu=({output.sig_ut}); for i in \"${{!arr[@]}}\"; do a=\"${{arr[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > \"${{orr[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orrd[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orru[$i]}}\"; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orr[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orru[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip >> \"${{orrd[$i]}}\"; else touch \"${{orr[$i]}}\" \"${{orrd[$i]}}\" \"${{orru[$i]}}\"; fi;done; for i in \"${{!arrt[@]}}\"; do a=\"${{arrt[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > \"${{orrt[$i]}}\"; cp \"${{orrt[$i]}}\" \"${{orrtd[$i]}}\"; cp \"${{orrt[$i]}}\" \"${{orrtu[$i]}}\"; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orrt[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orrtu[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[1]);if ($F[1] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip >> \"${{orrtd[$i]}}\"; else touch \"${{orrt[$i]}}\" \"${{orrtd[$i]}}\" \"${{orrtu[$i]}}\"; fi;done 2> {log}"
rule create_summary_snippet:
    input:  rules.run_DTU.output.res_t,
            rules.run_DTU.output.res_g,
            rules.run_DTU.output.res_p,
            rules.run_DTU.output.fig_F,
            rules.run_DTU.output.fig_P,
            rules.run_DTU.output.fig_PV,
            rules.run_DTU.output.fig_files,
            rules.themall.input.sig_g,
            rules.themall.input.sig_dg,
            rules.themall.input.sig_ug,
            rules.themall.input.sig_t,
            rules.themall.input.sig_dt,
            rules.themall.input.sig_ut,
            rules.themall.input.session
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DTU/{combo}/create_summary_snippet.log", combo=combo)
    conda:  "<REPO>/envs/"+DTUENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+DTUENV+""
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DTUENV} --loglevel DEBUG 2>> {log}"


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


