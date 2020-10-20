MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')
REFERENCE = subdict(config['MAPPING'],SETTINGS)['REFERENCE']
INDEX = str.split(',',subdict(config['MAPPING'],SETTINGS)['INDEX'])

rule generate_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/{unikey}/SAindex", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))
            tmp = temp(expand("TMP/{refdir}/ref.fa"),refdir=REFDIR),
            tmpa = temp(expand("TMP/{refdir}/ref.anno"),refdir=REFDIR)
    log:    expand("LOGS/{map}.idx.log", map=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'MAPPING')['ANNOTATION']]),
            genpath = expand("{refd}/{unikey}", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP')
    shell:  "rm -rf {params.tmpidx} && if [[ -f \"{params.genpath}\SAindex\" ]]; then ln -fs {params.genpath}\SAindex {output.idx} && touch {output.tmp} {output.tmpa} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {output.uidx} --outTmpDir {params.tmpidx} --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} 2> {log} && ln -s {params.genpath}/SAindex {output.idx} && cat {params.tmpidx}Log.out >> {log} && rm -f {params.tmpidx}Log.out && rm -rf {params.tmpidx};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{file}_unmapped_R1.fastq.gz",
                unmapped_r2 = "UNMAPPED/{file}_unmapped_R2.fastq.gz",
                tmp = temp("TMP/STAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = expand("{refd}/{unikey}", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
                #genpath = lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/{gen}{name}_{extension}".format(ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2)))
                anno = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION'],
                #anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && mv {output.tmp}.Unmapped.out.mate1 {output.unmapped_r1} 2>> {log} && mv {output.tmp}.Unmapped.out.mate2 {output.unmapped_r2} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                tmp = temp("TMP/STAROUT/{file}")
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                genpath = expand("{refd}/{unikey}", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
                #genpath = lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/{gen}{name}_{extension}".format(ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2)))
                anno = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION'],
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && mv {output.tmp}.Unmapped.out.mate* {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
