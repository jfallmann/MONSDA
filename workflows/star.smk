MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE),
    output: idx = expand("{ref}/{{dir}}/{map}/{{extension}}/{{gen}}{{name}}_{{extension}}_{map}.idx", ref=REFERENCE, map=MAPPERENV),
            tmp = temp("TMP/{dir}/{gen}{name}_{extension}.fa"),
            tmpa = temp("TMP/{dir}/{gen}{name}_{extension}.anno")
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{extension}}_{map}.idx.log", map=MAPPERENV)
    conda:  "snakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'MAPPING')['ANNOTATION']]),
            genpath = lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=wildcards.dir, map=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2)))+os.sep,
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
    shell:  "rm -rf {params.tmpidx} && if [[ -f \"{params.genpath}SAindex\" ]]; then ln -fs {params.genpath}SAindex {output.idx} && touch {output.tmp} {output.tmpa} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {params.tmpidx} --outTmpDir {params.tmpidx} --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} 2> {log} && ln -s {params.genpath}SAindex {output.idx} && cat {params.tmpidx}Log.out >> {log} && rm -f {params.tmpidx}Log.out && rm -rf {params.tmpidx};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped_r1 = "UNMAPPED/{file}_unmapped_r1.fastq.gz",
                unmapped_r2 = "UNMAPPED/{file}_unmapped_r2.fastq.gz",
                tmp = temp("TMP/STAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = lambda wildcards: "{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=source_from_sample(wildcards.file,config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && mv {output.tmp}.Unmapped.out.mate1 {output.unmapped_r1} 2>> {log} && mv {output.tmp}.Unmapped.out.mate2 {output.unmapped_r2} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                tmp = temp("TMP/STAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = lambda wildcards: "{ref}/{dir}/{map}/{extension}/".format(ref=REFERENCE, dir=source_from_sample(wildcards.file,config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING',2)),
                anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && mv {output.tmp}.Unmapped.out.mate* {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
