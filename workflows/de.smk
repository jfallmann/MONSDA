COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')
DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,'DE')

rule all:
    input:  expand("DE/{file}.csv",file=samplecond(SAMPLES,config))
            expand("COUNTS/{file}.summary", file=samplecond(SAMPLES,config)),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted.counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config)),
            expand("COUNTS/{file}_DONE",file=samplecond(SAMPLES,config))

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
    input: "SORTED_MAPPED/{file}_mapped_sorted.bam"
    output: "COUNTS/Featurecounter/{file}_mapped_sorted.counts",
            temp("COUNTS/Featurecounter/{file}.anno")
    log:    "LOGS/{file}/featurecount.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else ''
    shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} -a {output[1]} -o {output[0]} {input[0]} 2> {log}"

rule featurecount_unique:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts",
            temp("COUNTS/Featurecounter/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else ''
    shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} -a {output[1]} -o {output[0]} {input[0]} 2> {log}"

rule summarize_counts:
    input:  f = rules.count_fastq.output,
            m =rules.count_mappers.output,
            u =rules.count_unique_mappers.output
    log:    "LOGS/{file}/summarize_counts.log"
    output: "COUNTS/{file}.summary"
    conda:  "snakes/envs/base.yaml"
    threads: 1
    shell:  "arr=({input.f}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"

rule summarize_all:
    input:  c = expand(rules.summarize_counts.output, file=samplecond(SAMPLES,config)),
            f1 = expand(rules.featurecount.output, file=samplecond(SAMPLES,config)),
            f2 = expand(rules.featurecount_unique.output, file=samplecond(SAMPLES,config))
    output: a = "COUNTS/Features",
            u = "COUNTS/Features_unique",
            c = "COUNTS/Summary",
            t = expand("COUNTS/{file}_DONE",file=samplecond(SAMPLES,config))
    conda:  "snakes/envs/base.yaml"
    threads: 1
    params: bins = BINS
    shell:  "for i in {input.c};do if [[ $i == *\".summary\"*  ]];then cat $i >> {output.c};fi;done && for i in {input.f1};do if [[ $i == *\".counts\"*  ]];then cat $i\.summary >> {output.a};fi;done && for i in {input.f2};do if [[ $i == *\".counts\"*  ]];then cat $i\.summary >> {output.u};fi;done && touch {output.t}"

rule prepare_count_table:
    input:   cnt = expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/{file}.tbl"
    log:     "LOGS/DE/{file}.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  depara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'DE')['OPTIONS'][0].items()),
             bins = BINS
    shell: "python2 {params.bins}/DE/build_DESeq_table.py -l {input.cnt} -n > {output.tbl} 2> {log}"

rule run_deseq2:
    input:  cnt = expand(rules.prepare_count_table.output, file=samplecond(SAMPLES,config))
    output: csv = "DE/DESEQ2/{file}.csv"
    log:    "LOGS/DE/{file}.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params: debin = DEBIN,
            depara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'DE')['OPTIONS'][1].items()),
            bins = BINS
    shell: "Rscript {params.bins}/DE/ESeq2_diffexp.R {input.cnt} {output.csv} 2> {log}"

#rule themall:
#    input:  rules.summarize_counts.output
#    output: "COUNTS/DONE"
#    conda:  "snakes/envs/base.yaml"
#    threads: 1
#    params: bins = BINS
#    shell:  "touch {output}"

onsuccess:
    print("Workflow finished, no error")
