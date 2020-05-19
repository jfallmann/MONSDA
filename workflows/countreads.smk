COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

rule themall:
    input:  expand("COUNTS/Featurecounts_{feat}s/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config), feat=config['COUNTING']['FEATURES'].keys()),
            expand("COUNTS/Featurecounts_{feat}s/{file}_mapped_sorted.counts", file=samplecond(SAMPLES,config), feat=config['COUNTING']['FEATURES'].keys()),
            expand("COUNTS/{file}.summary", file=samplecond(SAMPLES,config))

if paired == 'paired':
    rule count_fastq:
        input:  r1 = lambda wildcards: expand("FASTQ/{rawfile}_{{read}}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),read=['R1','R2']),
                r2 = expand("TRIMMED_FASTQ/{{file}}_{read}_trimmed.fastq.gz", read=['R1','R2'])
        output: r1 = expand("COUNTS/{{file}}_raw_{read}_fq.count",read=['R1','R2']),
                r2 = expand("COUNTS/{{file}}_trimmed_{read}_fq.count",read=['R1','R2'])
        log:    "LOGS/{file}/countfastq.log"
        conda:  "snakes/envs/base.yaml"
        threads: 1
        shell:  "arr=({input.r1}); orr=({output.r1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log} && arr=({input.r2}); orr=({output.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log}"

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
    input:  m = "MAPPED/{file}_mapped_sorted.bam"
    output: m = "COUNTS/{file}_mapped.count"
    log:    "LOGS/{file}/countmappers.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T TMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_mappers:
    input:  u = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: u = "COUNTS/{file}_mapped_unique.count"
    log:    "LOGS/{file}/count_unique_mappers.log"
    conda:  "snakes/envs/samtools.yaml"
    threads: MAXTHREAD
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S 25% -T TMP -u |wc -l > {output.u} ;done 2>> {log}"

rule featurecount:
    input:  s = "MAPPED/{file}_mapped_sorted.bam",
    output: t = temp("COUNTS/Featurecounts_{feat}s/{file}_tmp.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{file}_mapped_sorted.counts"
    log:    "LOGS/{file}/featurecount_{feat}s.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.reads} 2> {log} && head -n2 {output.t} > {output.c} && tail -n+3 {output.t}|sort -k1,1 -k2,2n -k3,3n -u >> {output.c}"

rule featurecount_unique:
    input:  u = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
    output: t = temp("COUNTS/Featurecounts_{feat}s/{file}_tmp_uni.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecount_{feat}s_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.reads} 2> {log} && head -n2 {output.t} > {output.c} && tail -n+3 {output.t}|sort -k1,1 -k2,2n -k3,3n -u >> {output.c}"

rule summarize_counts:
    input:  f = rules.count_fastq.output,
            m = rules.count_mappers.output,
            u = rules.count_unique_mappers.output
    output: "COUNTS/{file}.summary"
    log:    "LOGS/{file}/summarize_counts.log"
    conda:  "snakes/envs/base.yaml"
    threads: 1
    shell:  "arr=({input.f}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"
