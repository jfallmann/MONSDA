include: "header.smk"

NAME=config["NAME"]["genomic"]

rule all:
    input:  expand("DONE/BED/{file}_{type}",file=samplecond(SAMPLES,config), type=['sorted','unique'])

rule LinkBeds:
    input:  "UCSC/{file}_mapped_{type}.bed"
    output: "BED/{file}_{type}.bed.gz"
    log:    "LOGS/Bed/linkbed{file}_{type}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed')
    shell:  "ln -s {params.abs} {output}"

rule AnnotateBed:
    input:  rules.LinkBeds.output
    output: "BED/{file}_anno_{type}.bed.gz"
    log:    "LOGS/Bed/annobeds_{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
            bins=BINS,
            anno=lambda wildcards: anno_from_file(wildcards.file, config),
            annop=config["ANNOTATE"]
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input[0]} -a {params.anno} {params.annop} |gzip > {output[0]}"

rule AddSequenceToBed:
    input:  rules.AnnotateBed.output
    output: "BED/{file}_seq_{type}.bed.gz",
            temp("BED/{file}_seq_{type}.tmp")
    log:    "LOGS/Beds/seq2bed_{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
    shell:  "fastaFromBed -fi {params.fasta} -bed {input[0]} -name+ -tab -s -fullHeader -fo {output[1]} && cut -d$'\t' -f2 {output[1]}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input[0]}) -|gzip  > {output[0]}"

rule MergeAnnoBed:
    input:  rules.AddSequenceToBed.output
    output: "BED/{file}_anno_{type}_merged.bed.gz"
    log:    "LOGS/Bed/mergebeds_{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=config["NAME"]),
            bins=BINS,
            anno=lambda wildcards: anno_from_file(wildcards.file, config)
    shell:  "zcat {input[0]}|perl -wlane 'print join(\"\t\",@F[0..5],$F[-2],$F[-1])' |bedtools merge -s -c 7,8 -o distinct -delim \"|\" |sort -k1,1 -k2,2n|gzip > {output[0]}"

rule themall:
    input:  rules.AnnotateBed.output
    output: "DONE/BED/{file}_{type}"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
