import glob, os, snakemake, inspect #, datetime

#snakemake -k --use-conda --ri --latency-wait 120 -j 16 -s Workflow/workflows/mapping.smk --configfile Workflow/config_dicty.json --directory ${PWD} -k
### optional with date
###snakemake --use-conda --ri --latency-wait 120 -j 16 -s mapping.smk --configfile config_dicty.json --directory ${PWD}/`date +%d-%m-%Y` -k

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Collection import *

QC=config["QC"]
ADAPTERS=config["ADAPTERS"]
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["NAME"]
BINS=config["BINS"]
SOURCE=sources(config)
SAMPLES=samples(config)

rule all:
    input:  "COUNTS/ENDS/Collect",
            "COUNTS/ENDS/CollectFQ",
            expand("COUNTS/{file}/Counts",file=SAMPLES)

rule count_fastq:
    input:  "FASTQ/{file}.fastq.gz",
            "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            "FILTERED/{file}_filtered.fastq.gz"
    output: "COUNTS/{file}_raw_fq.count",
            "COUNTS/{file}_trimmed_fq.count",
            "COUNTS/{file}_filtered_fq.count"
    conda:  "../envs/base.yaml"
    shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done"

rule count_mappers:
    input:  expand("SORTED_MAPPED/{{file}}_mapped_{cond}_sorted.bam", cond=config["COND4COUNT"]),
            expand("UNIQUE_MAPPED/{{file}}_mapped_{cond}_sorted_unique.bam", cond=config["COND4COUNT"]),
            expand("PHASED_MAPPED/{{file}}_mapped_{cond}_sorted_phased.bam", cond=config["COND4COUNT"])
    output: expand("COUNTS/{{file}}_mapped_{cond}.count", cond=config["COND4COUNT"]),
            expand("COUNTS/{{file}}_mapped_{cond}_unique.count", cond=config["COND4COUNT"]),
            expand("COUNTS/{{file}}_mapped_{cond}_phased.count", cond=config["COND4COUNT"])
    conda:  "../envs/samtools.yaml"
    threads: 20
    shell:  "export LC_ALL=C; arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > ${{orr[$i]}};done"

rule count_filtered:
    input:  "FILTERED/{file}_Gfiltered.sam.gz",
            "FILTERED/{file}_GPfiltered.sam.gz"
    output: "COUNTS/{file}_Gfiltered.count",
            "COUNTS/{file}_GPfiltered.count"
    conda:  "../envs/base.yaml"
    shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do zcat ${{arr[$i]}}| cut -d$'\t' -f1|sort -u | wc -l > ${{orr[$i]}};done"

rule count_ends:
    input:  expand("PHASED_MAPPED/{{file}}_mapped_{cond}_sorted_phased.bam", cond=config["COND4COUNT"]),
            expand("UNIQUE_MAPPED/{{file}}_mapped_{cond}_sorted_unique.bam", cond=config["COND4COUNT"]),
            expand("SORTED_MAPPED/{{file}}_mapped_{cond}_sorted.bam", cond=config["COND4COUNT"]),
    output: expand("COUNTS/ENDS/{{file}}_mapped_{cond}_sorted_phased.bam.ends.gz", cond=config["COND4COUNT"]),
            expand("COUNTS/ENDS/{{file}}_mapped_{cond}_sorted_unique.bam.ends.gz", cond=config["COND4COUNT"]),
            expand("COUNTS/ENDS/{{file}}_mapped_{cond}_sorted.bam.ends.gz", cond=config["COND4COUNT"])
    log:    "LOGS/CountEnds/{file}.log"
    conda: "../envs/base.yaml"
    threads: 1
    params: outdir = lambda w: "COUNTS/ENDS/"+str(os.path.dirname(w.file)),
            bins = BINS,
#            logfile= lambda w: "LOGS/CountEnds/"+str(w.file)+".log",
            offset = lambda wildcards: count_params(wildcards.file, config),
            clusterinfo = lambda wildcards: "{ref}/{gen}.tRNAscan_clusterInfo.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config))
    shell:  "arr=({input}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do bamfile=${{arr[$i]}}; logfile=${{bamfile##*/}}; python3 {params.bins}/Analysis/CountEnds.py -b ${{bamfile}} -o {params.outdir} -t {params.offset} -f {params.clusterinfo} -z {threads} &> LOGS/CountEnds/${{logfile}};done"

rule count_fastq_ends:
    input:  "FASTQ/{file}.fastq.gz",
            "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            "FILTERED/{file}_filtered.fastq.gz"
    output: "COUNTS/ENDS/{file}_raw_fq.ends.gz",
            "COUNTS/ENDS/{file}_trimmed_fq.ends.gz",
            "COUNTS/ENDS/{file}_filtered_fq.ends.gz"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bins = BINS,
    shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do perl {params.bins}/Analysis/CountFastqEnds.pl ${{arr[$i]}} | gzip > ${{orr[$i]}} ;done"

rule summarize_ends:
    input:  expand("COUNTS/ENDS/{file}_raw_fq.ends.gz", file=SAMPLES),
            expand("COUNTS/ENDS/{file}_trimmed_fq.ends.gz", file=SAMPLES),
            expand("COUNTS/ENDS/{file}_filtered_fq.ends.gz", file=SAMPLES)
    output: "COUNTS/ENDS/CollectFQ"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do {params.bins}/Shells/printFQEnds.sh ${{arr[$i]}} ${{orr}};done"

rule summarize_counts:
    input:  expand("COUNTS/{{file}}_mapped_{cond}.count", cond=config["COND4COUNT"]),
            expand("COUNTS/{{file}}_mapped_{cond}_unique.count", cond=config["COND4COUNT"]),
            expand("COUNTS/{{file}}_mapped_{cond}_phased.count", cond=config["COND4COUNT"]),
            expand("COUNTS/{{file}}_Gfiltered.count"),
            expand("COUNTS/{{file}}_GPfiltered.count"),
            expand("COUNTS/{{file}}_raw_fq.count"),
            expand("COUNTS/{{file}}_trimmed_fq.count"),
            expand("COUNTS/{{file}}_filtered_fq.count")
    output: "COUNTS/{file}/Counts"
    conda:  "../envs/base.yaml"
    threads: 1
    params: current = lambda w: w.file
    shell:  "arr=({input}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> COUNTS/{params.current}/Counts && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> COUNTS/{params.current}/Counts; else echo '0' >> COUNTS/{params.current}/Counts;fi;done"

rule themall:
    input:  expand("COUNTS/ENDS/{file}_mapped_cluster_sorted_phased.bam.ends.gz", file=SAMPLES),
            expand("COUNTS/ENDS/{file}_mapped_cluster_sorted_unique.bam.ends.gz", file=SAMPLES),
            expand("COUNTS/ENDS/{file}_mapped_cluster_sorted.bam.ends.gz", file=SAMPLES)
    output: "COUNTS/ENDS/Collect"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "arr=({input}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do {params.bins}/Shells/printEnds.sh ${{arr[$i]}} {output};done"
