include: "header.smk"

rule all:
    input:  expand("DONE/BED/{file}_{type}",file=samplecond(SAMPLES,config), type=['sorted','unique'])

checklist = list()
checklist2 = list()
for file in samplecond(SAMPLES,config):
    for type in ['sorted','unique']:
        checklist.append(os.path.isfile(os.path.abspath('UCSC/'+file+'_mapped_'+type+'.bed.gz')))
        checklist2.append(os.path.isfile(os.path.abspath('PEAKS/'+file+'_mapped_'+type+'.bed.gz')))

if all(checklist):
    rule BamToBed:
        input:  "UCSC/{file}_mapped_{type}.bed.gz"
        output: "BED/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('UCSC/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
elif all(checklist2):
    rule BamToBed:
        input:  "PEAKS/{file}_mapped_{type}.bed.gz"
        output: "BED/{file}_mapped_{type}.bed.gz"
        log:    "LOGS/Bed/linkbed{file}_{type}.log"
        conda:  "../envs/base.yaml"
        threads: 1
        params: abs = lambda wildcards: os.path.abspath('PEAKS/'+wildcards.file+'_mapped_'+wildcards.type+'.bed.gz')
        shell:  "ln -s {params.abs} {output}"
else:
    rule BamToBed:
        input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
                "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
        output: "BED/{file}_mapped_sorted.bed.gz",
                "BED/{file}_mapped_unique.bed.gz"
        log:    "LOGS/Bed/createbed{file}.log"
        conda:  "../envs/bedtools.yaml"
        threads: 1
        shell:  "bedtools bamtobed -i {input[0]} |gzip > {output[0]} && bedtools bamtobed -i {input[1]} |gzip > {output[1]}"

rule AnnotateBed:
    input:  rules.BamToBed.output
    output: "BED/{file}_anno_{type}.bed.gz"
    log:    "LOGS/Bed/annobeds_{type}_{file}.log"
    conda:  "../envs/perl.yaml"
    threads: 1
    params: bins=BINS,
            anno=lambda wildcards: anno_from_file(wildcards.file, config, 'annotation'),
            annop=config["ANNOTATE"],
            annof= lambda wildcards: "-s {feat}".format(feat=config["ANNOFEATURE"]) if config["ANNOFEATURE"] is not '' else ''
    shell:  "perl {params.bins}/Universal/AnnotateBed.pl -b {input[0]} -a {params.anno} {params.annof} {params.annop} |gzip > {output[0]}"

rule AddSequenceToBed:
    input:  rules.AnnotateBed.output
    output: "BED/{file}_anno_seq_{type}.bed.gz",
            temp("BED/{file}_anno_seq_{type}.tmp")
    log:    "LOGS/Bed/seq2bed_{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
    shell:  "export LC_ALL=C; fastaFromBed -fi {params.fasta} -bed <(zcat {input[0]}|cut -d$'\t' -f 1-6) -name+ -tab -s -fullHeader -fo {output[1]} && cut -d$'\t' -f2 {output[1]}|sed 's/t/u/ig'|paste -d$'\t' <(zcat {input[0]}) -|sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n |gzip  > {output[0]}"

rule MergeAnnoBed:
    input:  rules.AddSequenceToBed.output
    output: "BED/{file}_anno_seq_{type}_merged.bed.gz"
    log:    "LOGS/Bed/mergebeds_{type}_{file}.log"
    conda:  "../envs/bedtools.yaml"
    threads: 1
    params: fasta = lambda wildcards: "{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config)),
            bins=BINS,
            anno=lambda wildcards: anno_from_file(wildcards.file, config, 'annotation')
    shell:  "export LC_ALL=C; zcat {input[0]}|perl -wlane 'print join(\"\t\",@F[0..6],$F[-3],$F[-2])' |bedtools merge -s -c 7,8,9 -o distinct -delim \"|\" |sort --parallel={threads} -S 25% -T SORTTMP -t$'\t' -k1,1 -k2,2n|gzip > {output[0]}"

rule themall:
    input:  rules.MergeAnnoBed.output
    output: "DONE/BED/{file}_{type}"
    run:
        for f in output:
            with open(f, "w") as out:
                        out.write("DONE")
onsuccess:
    print("Workflow finished, no error")
