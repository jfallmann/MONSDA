COUNTBIN, COUNTENV = env_bin_from_config(config, 'COUNTING')

if not rundedup:
    rule themall:
        input:  expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),                
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/CPM_COUNTS_mapped_unique.counts.gz", file=samplecond(SAMPLES, config), feat='gene', combo=combo) if 'gene' in config['COUNTING']['FEATURES'].keys() else [],
                expand("COUNTS/Featurecounts_{feat}s/{combo}/TPM_COUNTS_mapped_unique.counts.gz", file=samplecond(SAMPLES, config), feat='gene', combo=combo) if 'gene' in config['COUNTING']['FEATURES'].keys() else [],
                expand("COUNTS/{combo}/{file}.summary", file=samplecond(SAMPLES, config), combo=combo)                

else:
    rule themall:
        input:  expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat=config['COUNTING']['FEATURES'].keys(), combo=combo),
                expand("COUNTS/Featurecounts_{feat}s/{combo}/CPM_COUNTS_mapped_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat='gene' if 'gene' in config['COUNTING']['FEATURES'].keys() else None, combo=combo) if 'gene' in config['COUNTING']['FEATURES'].keys() else [],
                expand("COUNTS/Featurecounts_{feat}s/{combo}/TPM_COUNTS_mapped_unique_dedup.counts.gz", file=samplecond(SAMPLES, config), feat='gene' if 'gene' in config['COUNTING']['FEATURES'].keys() else None, combo=combo) if 'gene' in config['COUNTING']['FEATURES'].keys() else [],
                expand("COUNTS/{combo}/{file}.summary", file=samplecond(SAMPLES, config), combo=combo)

if paired == 'paired':
    rule count_fastq:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
        output: r1 = "COUNTS/{combo}/{file}_R1_raw_fq.count",
                r2 = "COUNTS/{combo}/{file}_R2_raw_fq.count"
        log:    "LOGS/{combo}/{file}/countfastq.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        shell:  "arr=({input.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r1};done 2>> {log} && arr=({input.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r2};done 2>> {log}"

    rule count_trimmed_fastq:
        input:  r1 = expand("TRIMMED_FASTQ/{combo}/{{file}}_R1_trimmed.fastq.gz", combo=scombo),
                r2 = expand("TRIMMED_FASTQ/{combo}/{{file}}_R2_trimmed.fastq.gz", combo=scombo)
        output: r1 = "COUNTS/{combo}/{file}_R1_trimmed_fq.count",
                r2 = "COUNTS/{combo}/{file}_R2_trimmed_fq.count"
        log:    "LOGS/{combo}/{file}/countfastq.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        shell:  "arr=({input.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r1};done 2>> {log}; arr=({input.r2}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > {output.r2};done 2>> {log}"

else:
    rule count_fastq:
        input:  r1 = lambda wildcards: expand("FASTQ/{rawfile}.fastq.gz", rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: r1 = "COUNTS/{combo}/{file}_raw_fq.count"
        log:    "LOGS/{combo}/{file}/countfastq.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        shell:  "arr=({input.r1}); orr=({output.r1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log}"

    rule count_trimmed_fastq:
        input:  r1 = expand("TRIMMED_FASTQ/{combo}/{{file}}_trimmed.fastq.gz", combo=scombo)
        output: r1 = "COUNTS/{combo}/{file}_trimmed_fq.count"
        log:    "LOGS/{combo}/{file}/count_trimmedfastq.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        shell:  "arr=({input.r1}); orr=({output.r1}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do a=$(zcat ${{arr[$i]}}|wc -l ); echo $((a/4)) > ${{orr[$i]}};done 2>> {log}"

rule count_mappers:
    input:  m = expand("MAPPED/{combo}/{{file}}_mapped_sorted.bam", combo=scombo)
    output: m = "COUNTS/{combo}/{file}_mapped.count"
    log:    "LOGS/{combo}/{file}/countmappers.log"
    conda:  "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    params: sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}G -T TMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_mappers:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique.bam", combo=scombo)
    output: u = "COUNTS/{combo}/{file}_mapped_unique.count"
    log:    "LOGS/{combo}/{file}/count_unique_mappers.log"
    conda:  "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    params: sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}G -T TMP -u |wc -l > {output.u} ;done 2>> {log}"

rule count_dedup_mappers:
    input:  m = expand("MAPPED/{combo}/{{file}}_mapped_sorted_dedup.bam", combo=scombo)
    output: m = "COUNTS/{combo}/{file}_mapped_dedup.count"
    log:    "LOGS/{combo}/{file}/countdedupmappers.log"
    conda:  "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    params: sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "export LC_ALL=C; arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}G -T TMP -u |wc -l > {output.m} ;done 2>> {log}"

rule count_unique_dedup_mappers:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique_dedup.bam", combo=scombo)
    output: u = "COUNTS/{combo}/{file}_mapped_unique_dedup.count"
    log:    "LOGS/{combo}/{file}/count_unique_dedupmappers.log"
    conda:  "samtools.yaml"
    container: "oras://jfallmann/monsda:samtools"
    threads: MAXTHREAD
    params: sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "export LC_ALL=C; arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do samtools view -F 260 ${{arr[$i]}} | cut -d$'\t' -f1|sort --parallel={threads} -S {params.sortmem}G -T TMP -u |wc -l > {output.u} ;done 2>> {log}"

rule featurecount:
    input:  s = expand("MAPPED/{combo}/{{file}}_mapped_sorted.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_tmp.counts"),
            cts = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.s} 2> {log} && head -n2 {output.t} |gzip > {output.cts} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}G -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.cts} && mv {output.t}.summary {output.cts}.summary"

rule featurecount_unique:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_tmp_uni.counts"),
            cts = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_unique.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.u} 2> {log} && head -n2 {output.t} |gzip > {output.cts} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}G -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.cts} && mv {output.t}.summary {output.cts}.summary"

rule featurecount_dedup:
    input:  s = expand("MAPPED/{combo}/{{file}}_mapped_sorted_dedup.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_dedup_tmp.counts"),
            cts = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_dedup.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_dedup.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.s} 2> {log} && head -n2 {output.t} |gzip > {output.cts} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}G -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.cts} && mv {output.t}.summary {output.cts}.summary"

rule featurecount_unique_dedup:
    input:  u = expand("MAPPED/{combo}/{{file}}_mapped_sorted_unique_dedup.bam", combo=scombo)
    output: t = temp("COUNTS/Featurecounts_{feat}s/{combo}/{file}_dedup_tmp_uni.counts"),
            cts = "COUNTS/Featurecounts_{feat}s/{combo}/{file}_mapped_sorted_unique_dedup.counts.gz"
    log:    "LOGS/{combo}/{file}/featurecount_{feat}s_unique_dedup.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "COUNTING", COUNTENV)['OPTIONS'].get('COUNT', "")+' -t '+wildcards.feat+' -g '+config['COUNTING']['FEATURES'][wildcards.feat],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda w, resources: int(int(resources.mem_mb) / 1024)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.t} {input.u} 2> {log} && head -n2 {output.t} |gzip > {output.cts} && export LC_ALL=C; tail -n+3 {output.t}|sort --parallel={threads} -S {params.sortmem}G -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.cts} && mv {output.t}.summary {output.cts}.summary"

if rundedup:
    rule summarize_counts:
        input:  f = lambda wildcards: expand(rules.count_fastq.output, file=samplecond(SAMPLES, config), combo=combo) if paired == 'paired' else expand(rules.count_fastq.output, file=samplecond(SAMPLES, config), combo=combo),
                t = lambda wildcards: expand(rules.count_trimmed_fastq.output, file=samplecond(SAMPLES, config), combo=combo) if paired == 'paired' else expand(rules.count_trimmed_fastq.output, file=samplecond(SAMPLES, config), combo=combo),
                m = rules.count_mappers.output,
                u = rules.count_unique_mappers.output,
                d = rules.count_dedup_mappers.output,
                x = rules.count_unique_dedup_mappers.output
        output: "COUNTS/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/summarize_counts.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        params: curfile = lambda wildcards: wildcards.file
        shell:  "arr=({input.f}); filtered=();for val in ${{arr[@]}}; do [[ \"$val\" == *\"{params.curfile}\"* ]] && filtered+=(\"$val\");done; arr=filtered; alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" > {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.t}); alen=${{#arr[@]}}; filtered=();for val in ${{arr[@]}}; do [[ \"$val\" == *\"{params.curfile}\"* ]] && filtered+=(\"$val\");done; arr=filtered; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.d}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.x}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"

    rule prepare_count_table:
        input:  cts  = lambda wildcards: expand(rules.featurecount.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo),
                cts_u  = lambda wildcards: expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo),
                cts_d  = lambda wildcards: expand(rules.featurecount_dedup.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo),
                cts_ud  = lambda wildcards: expand(rules.featurecount_unique_dedup.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo)
        output:  tbl  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.counts.gz",
                 tbl_u  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.counts.gz",
                 tbl_d  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_dedup.counts.gz",
                 tbl_ud  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique_dedup.counts.gz",
                 anno  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.samples.gz",
                 anno_u = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.samples.gz",
                 anno_d = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_dedup.samples.gz",
                 anno_ud = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique_dedup.samples.gz"
        log:     "LOGS/DE/{combo}/Prepare_{feat}_count_table.log"
        conda:   "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        params:  samples = lambda wildcards, input: ','.join(input.cts),
                 samples_u = lambda wildcards, input: ','.join(input.cts_u),
                 samples_d = lambda wildcards, input: ','.join(input.cts_d),
                 samples_ud = lambda wildcards, input: ','.join(input.cts_ud),
                 bins = BINS
        shell:  "{params.bins}/Analysis/build_count_table_simple.py -r {params.samples} --table {output.tbl} --anno {output.anno} 2> {log} && {params.bins}/Analysis/build_count_table_simple.py -r {params.samples_u} --table {output.tbl_u} --anno {output.anno_u} 2> {log} && {params.bins}/Analysis/build_count_table_simple.py -r {params.samples_d} --table {output.tbl_d} --anno {output.anno_d} 2> {log} && {params.bins}/Analysis/build_count_table_simple.py -r {params.samples_ud} --table {output.tbl_ud} --anno {output.anno_ud} 2> {log}"

else:
    rule summarize_counts:
        input:  f = lambda wildcards: expand(rules.count_fastq.output, file=wildcards.file, combo=wildcards.combo), # samplecond(SAMPLES, config), combo=combo),
                t = lambda wildcards: expand(rules.count_trimmed_fastq.output,  file=wildcards.file, combo=wildcards.combo), # samplecond(SAMPLES, config), combo=combo),
                m = rules.count_mappers.output,
                u = rules.count_unique_mappers.output
        output: "COUNTS/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/summarize_counts.log"
        conda:  "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        params: curfile = lambda wildcards: wildcards.file
        shell:  "arr=({input.f}); filtered=();for val in ${{arr[@]}}; do [[ \"$val\" == *\"{params.curfile}\"* ]] && filtered+=(\"$val\");done; arr=filtered; alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" > {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.t}); filtered=();for val in ${{arr[@]}}; do [[ \"$val\" == *\"{params.curfile}\"* ]] && filtered+=(\"$val\");done; arr=filtered; alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.m}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done && arr=({input.u}); alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do echo -ne \"${{arr[$i]}}\t\" >> {output} && if [[ -s ${{arr[$i]}} ]]; then cat ${{arr[$i]}} >> {output}; else echo '0' >> {output};fi;done 2> {log}"

    rule prepare_count_table:
        input:  cts  = lambda wildcards: expand(rules.featurecount.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo),
                cts_u  = lambda wildcards: expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES, config), feat=wildcards.feat, combo=wildcards.combo)
        output:  tbl  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.counts.gz",
                 tbl_u  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.counts.gz",
                 anno  = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped.samples.gz",
                 anno_u = "COUNTS/Featurecounts_{feat}s/{combo}/COUNTS_mapped_unique.samples.gz"
        log:     "LOGS/DE/{combo}/Prepare_{feat}_count_table.log"
        conda:   "base.yaml"
        container: "oras://jfallmann/monsda:base"
        threads: 1
        params:  samples = lambda wildcards, input: ','.join(input.cts),
                 samples_u = lambda wildcards, input: ','.join(input.cts_u),                 
                 bins = BINS
        shell:  "{params.bins}/Analysis/build_count_table_simple.py -r {params.samples} --table {output.tbl} --anno {output.anno} 2> {log} && {params.bins}/Analysis/build_count_table_simple.py -r {params.samples_u} --table {output.tbl_u} --anno {output.anno_u} 2> {log}"
    
rule cpm_tpm_table:
    input:
        counts = rules.prepare_count_table.output.tbl_u if not rundedup else rules.prepare_count_table.output.tbl_ud
    output:
        cpm = "COUNTS/Featurecounts_{feat}s/{combo}/CPM_COUNTS_mapped_unique.counts.gz" if not rundedup else "COUNTS/Featurecounts_{feat}s/{combo}/CPM_COUNTS_mapped_unique_dedup.counts.gz",
        tpm = "COUNTS/Featurecounts_{feat}s/{combo}/TPM_COUNTS_mapped_unique.counts.gz" if not rundedup else "COUNTS/Featurecounts_{feat}s/{combo}/TPM_COUNTS_mapped_unique_dedup.counts.gz"
    log: "LOGS/DE/{combo}/CPM_TPM_{feat}_table.log"
    conda: "deseq2_DE.yaml"
    container: "oras://jfallmann/monsda:deseq2_DE"
    threads: 1
    params: bins = str.join(os.sep, [BINS, "Analysis/CPM_TPM_from_counts.R"]),
    shell: "Rscript {params.bins} {input.counts} {output.cpm} {output.tpm} 2> {log}"


