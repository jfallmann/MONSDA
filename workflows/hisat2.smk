MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = expand("{refd}/INDICES/{mape}_{unikey}/{pref}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]), pref=PREFIX),
            tmp = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            linkidx = lambda wildcards, output: str(os.path.abspath(str.join(os.sep,str(output.uidx[0]).split(os.sep)[:-1]))),
            tolink = lambda wildcards, output: str(os.path.abspath(str.join(os.sep,str(output.idx).split(os.sep)[:-1])))
    shell:  "if [[ -f \"{output.uidx}\" ]]; then ln -fs {params.linkidx} {output.idx} && touch {output.uidx} {output.tmp} && echo \"Found hisat index, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && {params.indexer}-build {params.ipara} -p {threads} {output.tmp} {output.uidx} 2> {log} && ln -fs {params.linkidx} {output.idx} && touch {output.uidx};fi"

hs2alg = MAPPERBIN.split(' ')[1] if ' ' in MAPPERBIN else MAPPERBIN

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                summary = "MAPPED/{file}.summary"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else '',
                pref = lambda wildcards, input: str.join(os.sep,[input.index,PREFIX]) if PREFIX != '' else input.index
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {params.pref} -1 {input.r1} -2 {input.r2} -S {output.mapped} --un-conc-gz {output.unmapped} --summary-file {output.summary} 2>> {log} && touch {output.unmapped}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
                summary = "MAPPED/{file}.summary"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else '',
                pref = lambda wildcards, input: str.join(os.sep,[input.index,PREFIX]) if PREFIX != '' else input.index
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {params.pref} -U {input.query} -S {output.mapped} --un-gz {output.unmapped} --summary-file {output.summary} 2>> {log} && touch {output.unmapped}"
