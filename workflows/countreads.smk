COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

#wildcard_constraints:
#    feat="!os.sep()"

log.warning('KEYS: '+str(list(config['COUNTING']['FEATURES'].keys())))

rule all:
    input:  #expand("COUNTS/Features_{feat}s", feat=config['COUNTING']['FEATURES'].keys()),
            #expand("COUNTS/Features_{feat}s_unique", feat=config['COUNTING']['FEATURES'].keys()),
            expand("COUNTS/{file}.summary", file=samplecond(SAMPLES,config)),
            "COUNTS/Summary",
            "COUNTS/DONE"

if paired == 'paired':
    rule count_fastq:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r3 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r4 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
        output: r1 = "COUNTS/{file}_raw_r1_fq.count",
                r2 = "COUNTS/{file}_raw_r2_fq.count",
                r3 = "COUNTS/{file}_trimmed_r1_fq.count",
                r4 = "COUNTS/{file}_trimmed_r2_fq.count"
        log:    "LOGS/{file}/countfastq.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        shell:  "arr=({input.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r1};done 2>> {log} && arr=({input.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r2};done 2>> {log} && arr=({input.r3}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r3};done 2>> {log} && arr=({input.r4}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r4};done 2>> {log}"

else:
    rule count_fastq:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        output: r1 = "COUNTS/{file}_raw_fq.count",
                r2 = "COUNTS/{file}_trimmed_fq.count"
        log:    "LOGS/{file}/countfastq.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        shell:  "arr=({input.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r1};done 2>> {log} && arr=({input.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r2};done 2>> {log}"

rule count_mappers:
    input:  m = "SORTED_MAPPED/{file}_mapped_sorted.bam"
    output: m = "COUNTS/{file}_mapped.count"
    log:    "LOGS/{file}/countmappers.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_mappers:
    input:  u = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: u = "COUNTS/{file}_mapped_unique.count"
    log:    "LOGS/{file}/count_unique_mappers.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T SORTTMP -u |wc -l > {output.u} ;done 2>> {log}"

rule featurecount:
    input:  s = "SORTED_MAPPED/{file}_mapped_sorted.bam",
    output: c = "COUNTS/Featurecounter_{feat}s/{file}_mapped_sorted.counts",
            t = temp("COUNTS/Featurecounter_{feat}s/{file}.anno")
    log:    "LOGS/{file}/featurecount_{feat}s.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output.t} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output.t} -o {output.c} {input.s} 2> {log}"

rule featurecount_unique:
    input:  u = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
    output: c = "COUNTS/Featurecounter_{feat}s/{file}_mapped_sorted_unique.counts",
            t = temp("COUNTS/Featurecounter_{feat}s/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_{feat}s_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output.t} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output.t} -o {output.c} {input.u} 2> {log}"

rule summarize_counts:
    input:  f = rules.count_fastq.output,
            m = rules.count_mappers.output,
            u = rules.count_unique_mappers.output
    output: "COUNTS/{file}.summary"
    log:    "LOGS/{file}/summarize_counts.log"
    conda:  "snakes/envs/base.yaml"
    threads: 1
    shell:  "arr=({input.f}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"

rule count_summary:
    input:  c = expand(rules.summarize_counts.output, file=samplecond(SAMPLES,config))
    output: c = "COUNTS/Summary"
    conda:  "snakes/envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input.c};do if [[ $i == *\".summary\"*  ]];then cat $i >> {output.c};fi;done"

rule themall:
    input:  f1 = expand(rules.featurecount.output.c, file=samplecond(SAMPLES,config),feat=config['COUNTING']['FEATURES'].keys()),
            f2 = expand(rules.featurecount_unique.output.c, file=samplecond(SAMPLES,config),feat=config['COUNTING']['FEATURES'].keys())
    output: a = "COUNTS/DONE"
    conda:  "snakes/envs/base.yaml"
    threads: 1
    params: bins = BINS,
            a = lambda x: expand("COUNTS/Features_{feat}s",feat=config['COUNTING']['FEATURES'].keys()),
            u = lambda x: expand("COUNTS/Features_{feat}s_unique",feat=config['COUNTING']['FEATURES'].keys())
    shell:  "for i in {input.f1};do if [[ $i == *\".counts\"*  ]];then cat $i\.summary >> {params.a};fi;done && for i in {input.f2};do if [[ $i == *\".counts\"*  ]];then cat $i\.summary >> {params.u};fi;done && touch {output.a}"

#f1 = expand("COUNTS/Featurecounter_{{feat}}/{file}_mapped_sorted.counts", file=samplecond(SAMPLES,config)),
            #f2 = expand("COUNTS/Featurecounter_{{feat}}/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config)),
            #,feat=list(config['COUNTING']['FEATURES'].keys()))#,feat=list(config['COUNTING']['FEATURES'].keys())),

onsuccess:
    print("Workflow finished, no error")

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
