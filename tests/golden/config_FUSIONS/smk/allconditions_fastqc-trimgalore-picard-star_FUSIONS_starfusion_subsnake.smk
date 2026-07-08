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

combo = 'fastqc-trimgalore-picard-star/starfusion'

scombo = 'fastqc-trimgalore-picard-star'

wildcard_constraints:
    combo = combo,
    scombo = scombo,
    read = "R1|R2",
    type = "sorted|sorted_unique" if not rundedup else "sorted|sorted_unique|sorted_dedup|sorted_unique_dedup"

FBIN, FENV = env_bin_from_config(config, 'FUSIONS')
if not 'star' in combo or not 'star' in scombo:
        log.warning('STAR-Fusion needs STAR chimeric output, can only be used with STAR in mapping step')
fusopts = tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS']
fastqmode = str(fusopts.get('FASTQ', "")).lower() in ("1", "true", "yes")
CTATLIB = fusopts.get('INDEX', "") or os.path.join(REFDIR, "CTAT", FENV)
if os.path.isfile(os.path.join(CTATLIB, "ref_annot.gtf")):
    CTATGENOMEDIR = CTATLIB
elif os.path.isfile(os.path.join(CTATLIB, "ctat_genome_lib_build_dir", "ref_annot.gtf")):
    CTATGENOMEDIR = os.path.join(CTATLIB, "ctat_genome_lib_build_dir")
else:
    CTATGENOMEDIR = None
if not rundedup:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}_fusions", combo=combo, file=samplecond(SAMPLES, config))
else:
    rule themall:
        input:  expand("FUSIONS/{combo}/{file}_{type}", combo=combo, file=samplecond(SAMPLES, config), type=['sorted', 'sorted_dedup'])
rule generate_ctat_lib:
    input:  fa = REFERENCE,
            anno = ANNOTATION
    output: lib = directory(CTATLIB),
            tmpfa = temp(expand("TMP/{fenv}/ctat_ref.fa", fenv=FENV)),
            tmpanno = temp(expand("TMP/{fenv}/ctat_ref.gtf", fenv=FENV))
    log:    expand("LOGS/{sets}/{fenv}.ctat.log", sets=SETS, fenv=FENV)
    conda:  "<REPO>/envs/"+FENV+".yaml"
    container: "oras://docker.io/jfallmann/monsda:"+FENV+""
    threads: MAXTHREAD
    params: sf = FBIN,
            bpara = lambda wildcards: tool_params(SAMPLES[0], None, config, "FUSIONS", FENV)['OPTIONS'].get('BUILD', "")
    shell:  "( zcat {input.fa} > {output.tmpfa} && zcat {input.anno} > {output.tmpanno} && mkdir -p {output.lib} && prep_genome_lib.pl --genome_fa {output.tmpfa} --gtf {output.tmpanno} --output_dir {output.lib} --CPU {threads} {params.bpara} ) &> {log}"
if CTATGENOMEDIR is not None:
    ctat_lib_input = CTATGENOMEDIR
else:
    ctat_lib_input = rules.generate_ctat_lib.output.lib
if not fastqmode:
    rule starfusion:
        input:  junction = expand("MAPPED/{scombo}/{{file}}.Chimeric.out.junction", scombo=scombo),
                lib = ctat_lib_input
        output: fusions = "FUSIONS/{combo}/{file}_fusions",
                normjunc = "FUSIONS/{combo}/{file}.Chimeric.norm.junction",
                outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
        log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
        conda:  "<REPO>/envs/"+FENV+".yaml"
        container: "oras://docker.io/jfallmann/monsda:"+FENV+""
        threads: MAXTHREAD
        params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                sf = FBIN
        shell:  "set +o pipefail; if [[ -s \"{input.junction}\" ]]; then ctat_chr=$(grep -v '^#' {input.lib}/ref_annot.gtf 2>/dev/null | head -1 | cut -f1 | grep -c '^chr' || true); junc_chr=$(awk '!/^#/ && $1!=\"chr_donorA\"{{print $1; exit}}' {input.junction} | grep -c '^chr' || true); if [[ \"$ctat_chr\" == \"1\" && \"$junc_chr\" == \"0\" ]]; then echo \"Adding chr prefix to junction to match CTAT lib\" >> {log}; awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print;next}} $1==\"chr_donorA\"{{print;next}} {{if($1!~/^chr/)$1=\"chr\"$1; if($4!~/^chr/)$4=\"chr\"$4; if($1==\"chrMT\")$1=\"chrM\"; if($4==\"chrMT\")$4=\"chrM\"; print}}' {input.junction} > {output.normjunc}; elif [[ \"$ctat_chr\" == \"0\" && \"$junc_chr\" == \"1\" ]]; then echo \"Stripping chr prefix from junction to match CTAT lib\" >> {log}; awk 'BEGIN{{OFS=\"\\t\"}} /^#/{{print;next}} $1==\"chr_donorA\"{{print;next}} {{sub(/^chr/,\"\",$1); sub(/^chr/,\"\",$4); if($1==\"M\")$1=\"MT\"; if($4==\"M\")$4=\"MT\"; print}}' {input.junction} > {output.normjunc}; else cp {input.junction} {output.normjunc}; fi; {params.sf} --genome_lib_dir {input.lib} -J {output.normjunc} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; else mkdir -p {output.outdir}; touch {output.normjunc}; echo \"File {input.junction} empty, no chimeric STAR output found\" >> {log}; fi; touch {output.fusions}"
else:
    if paired == 'paired':
        rule starfusion:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    lib = ctat_lib_input
            output: fusions = "FUSIONS/{combo}/{file}_fusions",
                    outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
            log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
            conda:  "<REPO>/envs/"+FENV+".yaml"
            container: "oras://docker.io/jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "{params.sf} --genome_lib_dir {input.lib} --left_fq {input.r1} --right_fq {input.r2} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; touch {output.fusions}"
    else:
        rule starfusion:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    lib = ctat_lib_input
            output: fusions = "FUSIONS/{combo}/{file}_fusions",
                    outdir = temp(directory("FUSIONS/{combo}/{file}_tmp"))
            log:    "LOGS/FUSIONS/{combo}/{file}_starfusion.log"
            conda:  "<REPO>/envs/"+FENV+".yaml"
            container: "oras://docker.io/jfallmann/monsda:"+FENV+""
            threads: MAXTHREAD
            params: fpara = lambda wildcards: tool_params(wildcards.file, None, config, "FUSIONS", FENV)['OPTIONS'].get('FUSION', ""),
                    sf = FBIN
            shell:  "{params.sf} --genome_lib_dir {input.lib} --left_fq {input.r1} --output_dir {output.outdir} --CPU {threads} {params.fpara} &> {log} && cp {output.outdir}/star-fusion.fusion_predictions.tsv {output.fusions} 2>> {log}; touch {output.fusions}"


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


