import glob, os, sys, inspect, snakemake

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
SOURCE=sources(config)
SAMPLES=samples(config)
if os.path.exists(SAMPLES[0]) is False:
    SAMPLES=sampleslong(config)

rule fastqc_raw:
    input: expand("FASTQ/{file}.fastq.gz",file=SAMPLES)
    output: report("QC/{file}_qc.zip", category="QC")
    wildcard_constraints:
        file="trimmed{0}"
    log:    "LOGS/{file}/fastqc_raw.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule fastqc_trimmed:
    input:  expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=SAMPLES),
            expand("QC/{file}_qc.zip", file=SAMPLES)
    output: report("QC/{file}_trimmed_qc.zip", category="QC")
    wildcard_constraints:
        file="trimmed{1}"
    log:    "LOGS/{file}/fastqc_trimmed.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input[0]} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule fastqc_mapped:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  report("QC/{file}_mapped_sorted_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule fastqc_uniquemapped:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"
