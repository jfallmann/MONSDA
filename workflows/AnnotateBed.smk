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
    shell:  "ln -s {input} {output}"

rule AddSequenceToBed:
    input:  rules.LinkBeds.output
    output: "BED/{file}_seq_{type}.bed.gz",
            temp("BED/{file}_seq_{type}.tmp")
    log:    "LOGS/Beds/seq2bed_{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=NAME),
    shell:  "fastaFromBed -fi {params.fasta} -bed {input[0]} -name+ -tab -s -fullHeader -fo {output[1]} && cut -d$'\t' -f2 {output[1]}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input[0]}) -|gzip  > {output[0]}"

rule AnnotateBed:
    input:  rules.AddSequenceToBed.output
    output: "BED/{file}_anno_{type}.bed.gz"
    log:    "LOGS/Bed/seq2beds_{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=config["NAME"]),
            bins=BINS,
            anno=lambda wildcards: anno_from_file(genome(wildcards.file, config), config)
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input} -a {params.anno} |gzip > {output}"

rule themall:
    input:  rules.AnnotateBed.output
    output: "DONE/BED/{file}_{type}"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
