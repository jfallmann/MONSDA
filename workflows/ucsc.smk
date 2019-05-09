import glob, os, sys, inspect, snakemake

###snakemake -n -j 20 --use-conda -s Workflow/workflows/mapping_paired.smk
###--configfile Workflow/config_compare.json --directory ${PWD}
###--printshellcmds 2> run.log

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Collection import *

ADAPTERS=config["ADAPTERS"]
REFERENCE=config["REFERENCE"]
GENOME=config["GENOME"]
NAME=config["UCSCNAME"]
BINS=config["BINS"]
SOURCE=sources(config)
SAMPLES=samples(config)

rule all:
    input:  expand("DONE/{file}{placeholder}tracks",file=SAMPLES, placeholder=get_placeholder(config))

rule bamtobed:
    input:  "SORTED_MAPPED/{file}_mapped{placeholder}sorted.bam",
            "UNIQUE_MAPPED/{file}_mapped{placeholder}sorted_unique.bam"
    output: "UCSC/{file}_mapped{placeholder}sorted.bed",
            "UCSC/{file}_mapped{placeholder}unique.bed"
    conda:  "../envs/bedtools.yaml"
    shell:  "bedtools bamtobed -i {input[0]} > {output[0]} && bedtools bamtobed -i {input[1]} > {output[1]} "

rule index_fa:
    input:  expand("{ref}/{{org}}/{{gen}}.{name}.fa",ref=REFERENCE, name=NAME),#"{ref}/{{org}}/{{gen}}_{name}.fa", ref=config["REFERENCE"], name =config["NAME"])#, gen=pathstogenomes(SAMPLES, config))
    output: expand("{ref}/{{org}}/{{gen}}.{name}.fa.fai",ref=REFERENCE, name=NAME)#"{ref}/{{org}}/{{gen}}_{name}.fa.fai", ref=config["REFERENCE"], name =config["NAME"])#, gen=pathstogenomes(SAMPLES, config))
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "for i in {input};do {params.bins}/Preprocessing/indexfa.sh $i 2> {log};done"

rule get_chromsize_genomic:
    input: expand("{ref}/{{org}}/{{gen}}.{name}.fa.fai",ref=REFERENCE, name=NAME)
        #expand("{ref}/{{org}}/{{gen}}_{name}.fa.fai", ref=config["REFERENCE"], name =config["NAME"])
    output: expand("{ref}/{{org}}/{{gen}}.{name}.chrom.sizes",ref=REFERENCE, name=NAME)
        #expand("{ref}/{{org}}/{{gen}}_{name}_chrom.sizes", ref=config["REFERENCE"], name =config["NAME"])
    conda:  "../envs/samtools.yaml"
    params: bins = BINS
    threads: 1
    shell:  "cut -f1,2 {input} |sed 's/^chr//g' > {output} 2> {log}"

rule BedToBedg:
    input:  "UCSC/{file}_mapped{placeholder}sorted.bed",
            "UCSC/{file}_mapped{placeholder}unique.bed"
    output: "UCSC/{file}_mapped{placeholder}sorted.fw.bedg.gz",
            "UCSC/{file}_mapped{placeholder}sorted.re.bedg.gz",
            "UCSC/{file}_mapped{placeholder}unique.fw.bedg.gz",
            "UCSC/{file}_mapped{placeholder}unique.re.bedg.gz"
    conda:  "../envs/perl.yaml"
    params: out=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
            bins=BINS,
            sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=NAME)
    shell:  "perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[0]} -c {params.sizes} -v on -x {output[0]} -y {output[1]} -a track && perl {params.bins}/Universal/Bed2Bedgraph.pl -f {input[1]} -c {params.sizes} -v on -x {output[2]} -y {output[3]} -a track"

### This step generates bigwig files for peaks which can then be copied to a web-browsable directory and uploaded to UCSC via the track field
rule BedgToUCSC:
    input:  "UCSC/{file}_mapped{placeholder}sorted.fw.bedg.gz",
            "UCSC/{file}_mapped{placeholder}sorted.re.bedg.gz",
            "UCSC/{file}_mapped{placeholder}unique.fw.bedg.gz",
            "UCSC/{file}_mapped{placeholder}unique.re.bedg.gz"
    output: "UCSC/{file}_mapped{placeholder}sorted.fw.bw",
            "UCSC/{file}_mapped{placeholder}sorted.re.bw",
            "UCSC/{file}_mapped{placeholder}unique.fw.bw",
            "UCSC/{file}_mapped{placeholder}unique.re.bw",
            temp("UCSC/{file}_{placeholder}_fw_tmp"),
            temp("UCSC/{file}_{placeholder}_re_tmp"),
            temp("UCSC/{file}_{placeholder}_unique_fw_tmp"),
            temp("UCSC/{file}_{placeholder}_unique_re_tmp")
    conda:  "../envs/ucsc.yaml"
    params: sizes = lambda wildcards: "{ref}/{gen}.{name}.chrom.sizes".format(ref=REFERENCE,gen=genomepath(wildcards.file,config),name=NAME)
    shell:  "export LC_ALL=C; zcat {input[0]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[4]} && bedGraphToBigWig {output[4]} {params.sizes} {output[0]} && zcat {input[1]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[5]} && bedGraphToBigWig {output[5]} {params.sizes} {output[1]} && zcat {input[2]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[6]} && bedGraphToBigWig {output[6]} {params.sizes} {output[2]} && zcat {input[3]} |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n > {output[7]} && bedGraphToBigWig {output[7]} {params.sizes} {output[3]}"

rule themall:
    input:  expand("UCSC/{{file}}_mapped{placeholder}sorted.fw.bw",placeholder=get_placeholder(config)),
            expand("UCSC/{{file}}_mapped{placeholder}sorted.re.bw",placeholder=get_placeholder(config)),
            expand("UCSC/{{file}}_mapped{placeholder}sorted.fw.bedg.gz",placeholder=get_placeholder(config)),
            expand("UCSC/{{file}}_mapped{placeholder}sorted.re.bedg.gz",placeholder=get_placeholder(config))
    output: expand("DONE/{{file}}{placeholder}tracks",placeholder=get_placeholder(config))
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
