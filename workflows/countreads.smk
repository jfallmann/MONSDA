include: "header.smk"

rule all:
    input:  expand("COUNTS/{file}/Counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/{rawfile}/Counts",rawfile=SAMPLES),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted.counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config)),
            "COUNTS/DONE"

if config['MAPPING'] is 'paired':
    rule count_fastq:
        input:  expand("FASTQ/{rawfile}_r1.fastq.gz", rawfile=SAMPLES),
                expand("FASTQ/{rawfile}_r2.fastq.gz", rawfile=SAMPLES),
                expand("TRIMMED_FASTQ/{rawfile}_r1_trimmed.fastq.gz", rawfile=SAMPLES),
                expand("TRIMMED_FASTQ/{rawfile}_r2_trimmed.fastq.gz", rawfile=SAMPLES)
        output: "COUNTS/{file}_raw_r1_fq.count",
                "COUNTS/{file}_raw_r2_fq.count",
                "COUNTS/{file}_trimmed_r1_fq.count",
                "COUNTS/{file}_trimmed_r2_fq.count"
        log:    "LOGS/{file}/countfastq.log"
        conda:  "../envs/base.yaml"
        threads: 1
        shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}} 2>> {log};done"

else:
    rule count_fastq:
        input:  expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES),
                expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=samplecond(SAMPLES,config))
        output: "COUNTS/{file}_raw_fq.count",
                "COUNTS/{file}_trimmed_fq.count"
        log:    "LOGS/{file}/countfastq.log"
        conda:  "../envs/base.yaml"
        threads: 1
        shell:  "arr=({input}); alen=${{#arr[@]}}; orr=({output}); for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}} 2>> {log} ;done"

rule count_mappers:
    input:  m = expand("SORTED_MAPPED/{file}_mapped_sorted.bam", file=samplecond(SAMPLES,config))
    output: m = "COUNTS/{file}_mapped.count"
    log:    "LOGS/{file}/countmappers.log"
    conda:  "../envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_mappers:
    input:  u = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: u = "COUNTS/{file}_mapped_unique.count"
    log:    "LOGS/{file}/countmappers.log"
    conda:  "../envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}};  for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > {output.u};done 2>> {log}"

rule featurecount:
    input:  expand("SORTED_MAPPED/{file}_mapped_sorted.bam", file=samplecond(SAMPLES,config))
    output: "COUNTS/Featurecounter/{file}_mapped_sorted.counts",
            temp("COUNTS/Featurecounter/{file}.anno")
    log:    "LOGS/{file}/featurecount.log"
    conda:  "../envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(REFERENCE,source_from_sample(wildcards.file).split(os.sep)[0],config["ANNOTATION"][source_from_sample(wildcards.file).split(os.sep)[0]]['counting'])),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")[0].items()),
            paired = lambda: '-p' if config['MAPPINGTYPE'] == 'paired' else ''
    shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} -a {output[1]} -o {output[0]} {input[0]} 2> {log}"

rule featurecount_unique:
    input:  expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: "COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts",
            temp("COUNTS/Featurecounterunique/{file}.anno")
    log:    "LOGS/{file}/featurecount_unique.log"
    conda:  "../envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(REFERENCE,source_from_sample(wildcards.file).split(os.sep)[0],config["ANNOTATION"][source_from_sample(wildcards.file).split(os.sep)[0]]['counting'])),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")[0].items()),
            paired = lambda: '-p' if config['MAPPINGTYPE'] == 'paired' else ''
    shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} -a {output[1]} -o {output[0]} {input[0]} 2> {log}"

rule summarize_counts:
    input:  rules.count_fastq.output,
            rules.count_mappers.output
    log:    "LOGS/{file}/summarize_counts.log"
    output: "COUNTS/{file}/Counts"
    conda:  "../envs/base.yaml"
    threads: 1
    params: current = lambda w,input: os.path.dirname(input[0])
    shell:  "arr=({input}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> COUNTS/{params.current}/Counts && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> COUNTS/{params.current}/Counts; else echo '0' >> COUNTS/{params.current}/Counts;fi;done 2 > {log}"

onsuccess:
    print("Workflow finished, no error")

rule themall:
    input:  expand(rules.summarize_counts.output, file=samplecond(SAMPLES,config)),
            expand(rules.featurecount.output, file=samplecond(SAMPLES,config)),
            expand(rules.featurecount_unique.output, file=samplecond(SAMPLES,config))
    output: "COUNTS/DONE"
    conda:  "../envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "touch {output}"

###rnacounter and cufflinks are to be added later
#rule RNAcountReads:
#   input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
#           "COUNTS/Featurecounter/{file}_mapped_sorted.counts"
#   output: "COUNTS/RNAcounter/{file}_mapped_sorted.counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.gene_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.transcript_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.exon_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.intron_counts"
#   shell:  "/usr/bin/rnacounter --nh -n 1 -t genes,transcripts,exons,introns {input[0]} {ANNOTATION} > {output[0]} && grep 'gene' {output[0]} > {output[1]} && grep 'transcript' {output[0]} > {output[2]} && grep 'exon' {output[0]} > {output[3]} && grep 'intron' {output[0]} > {output[4]}"
#
#rule RNAcountReads_uniq:
#   input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
#           "COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts"
#   output: "COUNTS/RNAcounter/{file}_mapped_sorted_unique.counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.gene_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.transcript_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.exon_counts",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.intron_counts"
#   shell:  "/usr/bin/rnacounter -n 1 -t genes,transcripts,exons,introns {input[0]} {ANNOTATION} > {output[0]} && grep 'gene' {output[0]} > {output[1]} && grep 'transcript' {output[0]} > {output[2]} && grep 'exon' {output[0]} > {output[3]} && grep 'intron' {output[0]} > {output[4]}"
#
#rule cufflinks:
#   input:  "SORTED_MAPPED/{file}_mapped_sorted.bam",
#           "COUNTS/RNAcounter/{file}_mapped_sorted.counts"
#   output: "QUANT/Cufflinks/{file}/transcripts.gtf"
#   log:    "LOGS/Cufflinks/{file}.log"
#   params: out="QUANT/Cufflinks/{file}"
#   threads: 20
#   shell:  "cufflinks -o {params.out} -p {threads} -G {ANNOTATION} {input[0]}"
#
#rule cufflinks_uniq:
#   input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
#           "COUNTS/RNAcounter/{file}_mapped_sorted_unique.counts"
#   output: "QUANT/Cufflinks/{file}_unique/transcripts.gtf"
#   log:    "LOGS/Cufflinks/{file}.log"
#   params: out="QUANT/Cufflinks/{file}_unique"
#   threads: 20
#   shell:  "cufflinks -o {params.out} -p {threads} -G {ANNOTATION} {input[0]}"
