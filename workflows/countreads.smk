COUNTBIN, COUNTENV = env_bin_from_config3(config, 'COUNTING')

if not rundedup:
    rule themall:
        input:  expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/{combo}/{file}.summary", file=samplecond(SAMPLES, config), combo=combo)
else:
    rule themall:
        input:  expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/{combo}/{file}.summary", file=samplecond(SAMPLES, config), combo=combo)

if paired == 'paired':
    rule count_fastq:
        input:  r1 = lambda wildcards: expand("FASTQ/{rawfile}_{{read}}.fastq.gz", rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = expand("TRIMMED_FASTQ/{combo}/{{file}}_{{read}}_trimmed.fastq.gz", combo=scombo)
        output: r1 = "COUNTS/{combo}/{file}_raw_{read}_fq.count",
                r2 = "COUNTS/{combo}/{file}_trimmed_{read}_fq.count"
        log:    "LOGS/{combo}/{file}_{read}/countfastq.log"
        conda:  "base.yaml"
        threads: 1
        shell:  "arr=({input.r1}); orr=({output.r1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log} && arr=({input.r2}); orr=({output.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log}"

else:
    rule count_fastq:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = expand("TRIMMED_FASTQ/{combo}/{{file}}_trimmed.fastq.gz", combo=scombo)
        output: r1 = "COUNTS/{combo}/{file}_raw_fq.count",
                r2 = "COUNTS/{combo}/{file}_trimmed_fq.count"
        log:    "LOGS/{combo}/{file}/countfastq.log"
        conda:  "base.yaml"
        threads: 1
        shell:  "arr=({input.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r1};done 2>> {log} && arr=({input.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r2};done 2>> {log}"

rule count_mappers:
    input:  m = expand("MAPPED/{combo}/{{file}}_mapped_sorted.bam", combo=scombo)
    output: m = "COUNTS/{combo}/{file}_mapped.count"
    log:    "LOGS/{combo}/{file}/countmappers.log"
    conda:  "samtools.yaml"
    threads: MAXTHREAD
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}% -T TMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_mappers:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique.bam", combo=scombo)
    output: u = "COUNTS/{combo}/{file}_mapped_unique.count"
    log:    "LOGS/{combo}/{file}/count_unique_mappers.log"
    conda:  "samtools.yaml"
    threads: MAXTHREAD
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}% -T TMP -u |wc -l > {output.u} ;done 2>> {log}"

rule count_dedup_mappers:
    input:  m = expand("MAPPED/{combo}/{{file}}_mapped_sorted_dedup.bam", combo=scombo)
    output: m = "COUNTS/{combo}/{file}_mapped_dedup.count"
    log:    "LOGS/{combo}/{file}/countdedupmappers.log"
    conda:  "samtools.yaml"
    threads: MAXTHREAD
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}% -T TMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_dedup_mappers:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique_dedup.bam", combo=scombo)
    output: u = "COUNTS/{combo}/{file}_mapped_unique_dedup.count"
    log:    "LOGS/{combo}/{file}/count_unique_dedupmappers.log"
    conda:  "samtools.yaml"
    threads: MAXTHREAD
    params: sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}% -T TMP -u |wc -l > {output.u} ;done 2>> {log}"

rule featurecount:
    input:  s = expand("MAPPED/{combo}/{{file}}_mapped_sorted.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_tmp.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.s} 2> {log} && head -n2 {output.t} |gzip > {output.c} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.c} && mv {output.t}.summary {output.c}.summary"

rule featurecount_unique:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_tmp_uni.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_unique.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.u} 2> {log} && head -n2 {output.t} |gzip > {output.c} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.c} && mv {output.t}.summary {output.c}.summary"

rule featurecount_dedup:
    input:  s = expand("MAPPED/{combo}/{{file}}_mapped_sorted_dedup.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_dedup_tmp.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_dedup.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_dedup.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.s} 2> {log} && head -n2 {output.t} |gzip > {output.c} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.c} && mv {output.t}.summary {output.c}.summary"

rule featurecount_unique_dedup:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique_dedup.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_dedup_tmp_uni.counts"),
            c = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique_dedup.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_unique_dedup.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.u} 2> {log} && head -n2 {output.t} |gzip > {output.c} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.c} && mv {output.t}.summary {output.c}.summary"

if rundedup:
    rule summarize_counts:
        input:  f = lambda wildcards: expand(rules.count_fastq.output, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0], file=samplecond(SAMPLES, config), read=["R1", "R2"], combo=combo) if paired == 'paired' else expand(rules.count_fastq.output, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0], file=samplecond(SAMPLES, config), combo=combo),
                m = rules.count_mappers.output,
                u = rules.count_unique_mappers.output,
                d = rules.count_dedup_mappers.output,
                x = rules.count_unique_dedup_mappers.output
        output: "COUNTS/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/summarize_counts.log"
        conda:  "base.yaml"
        threads: 1
        shell:  "arr=({input.f}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.d}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.x}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"
else:
    rule summarize_counts:
        input:  f = lambda wildcards: expand(rules.count_fastq.output, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0], file=samplecond(SAMPLES, config), read=["R1", "R2"], combo=combo) if paired == 'paired' else expand(rules.count_fastq.output, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0], file=samplecond(SAMPLES, config), combo=combo),
                m = rules.count_mappers.output,
                u = rules.count_unique_mappers.output
        output: "COUNTS/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/summarize_counts.log"
        conda:  "base.yaml"
        threads: 1
        shell:  "arr=({input.f}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"
