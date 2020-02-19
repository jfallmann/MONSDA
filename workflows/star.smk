MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE),
    output: idx = expand("{ref}/{{dir}}/{map}/{{extension}}/{{gen}}{{name}}_{{extension}}_{map}.idx", ref=REFERENCE, map=MAPPERENV),
            tmp = temp("TMPIDX/{dir}/{gen}{name}_{extension}.fa"),
            tmpa = temp("TMPIDX/{dir}/{gen}{name}_{extension}.anno")
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{extension}}_{map}.idx.log", map=MAPPERENV)
    conda:  "snakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'MAPPING')['ANNOTATION']]),
            genpath = lambda wildcards: "{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=wildcards.dir, map=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2))
    shell:  "if [[ -f \"{params.genpath}SAindex\" ]]; then touch {output.idx} {output.tmp} {output.tmpa} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && rm -rf TMPSTAR && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outTmpDir TMPSTAR --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} 2> {log} && touch {output.idx} && cat Log.out >> {log} && rm -f Log.out && rm -rf TMPSTAR ;fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                tmp = temp("TMPSTAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = lambda wildcards: "{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp} --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}Aligned.out.sam {output.mapped} 2>> {log} && cat  {output.tmp}Unmapped.out.mate* > {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp} 2>>{log}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                tmp = temp("TMPSTAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = lambda wildcards: "{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp} --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}Aligned.out.sam {output.mapped} 2>> {log} && cat  {output.tmp}Unmapped.out.mate* > {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp} 2>>{log}"
