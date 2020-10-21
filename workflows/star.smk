MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}/{unikey}/star.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
            tmp = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV)),
            tmpa = temp(expand("TMP/{mape}/ref.anno", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            anno = ANNO,
            linkidx = lambda wildcards, output: str(os.path.abspath(str.join(os.sep,[str(output.uidx[0]).split(os.sep)][:-1]))),
            genpath = expand("{refd}/INDICES/{mape}/{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP')
    shell:  "rm -rf {params.tmpidx} && if [[ -f \"{params.genpath}\SAindex\" ]]; then ln -fs {params.genpath} {output.idx} && touch {output.tmp} {output.tmpa} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {params.tmpidx} --outTmpDir {params.tmpidx} --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} 2> {log} && ln -s {params.linkidx} {output.idx} && touch output.uidx && cat {params.tmpidx}Log.out >> {log} && rm -f {params.tmpidx}Log.out && rm -rf {params.tmpidx};fi"

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
                genpath = expand("{refd}/INDICES/{mape}/{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
                anno = ANNO,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
                #genpath = lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/{gen}{name}_{extension}".format(ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), mape=MAPPERENV, extension=check_tool_params(SAMPLES[0], None ,config, 'MAPPING',2)))
                #anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'MAPPING')['ANNOTATION']]),
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
                genpath = expand("{refd}/INDICES/{mape}/{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0])),
                anno = ANNO,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && mv {output.tmp}.Unmapped.out.mate* {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
