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
NAME=config["NAME"] if config["NAME"] else ''
BINS=config["BINS"]
MAPPERENV=config["MAPPERENV"]
MAPPERBIN=config["MAPPERBIN"]
SOURCE=sources(config)
SAMPLES=samples(config)
if os.path.exists(SAMPLES[0]) is False:
    SAMPLES=sampleslong(config)

#include: "porechop.smk" #We lack a proper adapter trimming tool
include: "fastqc.smk"           #Need alternative for really long reads

rule all:
    input: expand("DONE/{file}_mapped",file=SAMPLES),
           "QC/Multi/DONE"

rule mapping:
    input:  "TRIMMED_FASTQ/{file}_trimmed.fastq.gz", "QC/{file}_trimmed_qc.zip" if "ON" in config["QC"] else "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: "{mappar}".format(mappar=mapping_params(wildcards.file, None ,config)),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME['genomic'])),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} {params.ref} {input[0]} | tee -a >(grep -v -P '\t4\t' >{output[0]}) | grep -P '\t4\t' > {output[1]} 2> {log}"

rule gzipsam:
    input: "MAPPED/{file}_mapped.sam",
           "UNMAPPED/{file}_unmapped.fastq",
    output: report("MAPPED/{file}_mapped.sam.gz", category="ZIPIT"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/gzipsam.log"
    conda:  "../envs/base.yaml"
    threads: 20
    shell: "pigz -k -p {threads} -f {input} > {output} 2> {log}"

rule sortsam:
    input:  "MAPPED/{file}_mapped.sam.gz"
    output: report("SORTED_MAPPED/{file}_mapped_sorted.sam.gz", category="SORTING"),
            temp("SORTED_MAPPED/{file}_mapped_header.gz"),
            temp("SORTTMP/{file}")
    log:    "LOGS/{file}/sortsam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    shell: "set +o pipefail;samtools view -H {input[0]}|grep -P '^@HD' |pigz -p {threads} -f > {output[1]} ; samtools view -H {input[0]}|grep -P '^@SQ'|sort -t$'\t' -k1,1 -k2,2V |pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@RG'|pigz -p {threads} -f >> {output[1]} ; samtools view -H {input[0]}|grep -P '^@PG'|pigz -p {threads} -f >> {output[1]} ; export LC_ALL=C;zcat {input[0]} | grep -v \"^@\"|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k3,3V -k4,4n - |pigz -p {threads} -f > {output[2]} ; cat {output[1]} {output[2]} > {output[0]} 2> {log}"

rule sam2bam:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz", "QC/{file}_mapped_sorted_qc.zip" if  "ON" in config["QC"] else "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output: report("SORTED_MAPPED/{file}_mapped_sorted.bam", category="2BAM"),
            "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    log:    "LOGS/{file}/sam2bam.log"
    conda: "../envs/samtools.yaml"
    threads: 20
    params: bins = BINS,
#            sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file)))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
        #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule uniqsam:
    input: "SORTED_MAPPED/{file}_mapped_sorted.sam.gz",
           "SORTED_MAPPED/{file}_mapped_sorted.bam",
           "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output: report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz", category="UNIQUE")
    log: "LOGS/{file}/uniqsam.log"
    conda: "../envs/base.yaml"
    threads: 20
    params: bins=BINS
    shell:  "{params.bins}/Shells/UniqueSam_woPicard.sh {input[0]} {output[0]} {threads} 2> {log}"

rule sam2bamuniq:
    input:   "UNIQUE_MAPPED/{file}_mapped_sorted_unique.sam.gz",
             "SORTED_MAPPED/{file}_mapped_sorted.bam",
             "SORTED_MAPPED/{file}_mapped_sorted.bam.bai"
    output:  report("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", category="2BAM"),
             "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    log:     "LOGS/{file}/sam2bamuniq.log"
    conda:   "../envs/samtools.yaml"
    threads: 20
    params: bins=BINS,
            #           sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE, gen=genomepath(wildcards.file,config), name=NAME),
            fn = lambda wildcards: "{fn}".format(fn=sample_from_path(wildcards.file))
    shell: "zcat {input[0]} | samtools view -bS - | samtools sort -T {params.fn} -o {output[0]} --threads {threads} && samtools index {output[0]} 2> {log}"
           #"{params.bins}/Shells/Sam2Bam.sh {input[0]} {params.sizes} {output}"

rule multiqc:
    input:  expand("QC/{file}_qc.zip", file=SAMPLES),
            expand("QC/{file}_trimmed_qc.zip", file=SAMPLES),
            expand("QC/{file}_mapped_sorted_qc.zip", file=SAMPLES),
            expand("QC/{file}_mapped_sorted_unique_qc.zip", file=SAMPLES)
    output: report("QC/Multi/DONE", category="QC")
    log:    "LOGS/multiqc.log"
    conda:  "../envs/qc.yaml"
    shell:  "OUT=$(dirname {output}); multiqc -k json -z -o $OUT $PWD 2> {log} && touch QC/Multi/DONE"

rule themall:
    input:  "QC/Multi/DONE", "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam" if "ON" in config["QC"] else "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "DONE/{file}_mapped"
    run:
        for f in output:
            with open(f, "w") as out:
                out.write("DONE")
