MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')
REFERENCE = subdict(config['MAPPING'],SETTINGS)['REFERENCE']
INDEX = str.split(',',subdict(config['MAPPING'],SETTINGS)['INDEX'])

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx1 = INDEX[0],
	        idx2 = INDEX[1],
            uidx1 = expand("{refd}/{unikey}.idx", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))
            uidx2 = expand("{refd}/{unikey}.idx", ref=REFDIR, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))+'second'
    log:    expand("LOGS/{map}.idx.log", map=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items())
    shell: "{params.indexer} --threads {threads} {params.ipara} -d {input.fa} -x {output.uidx1} -y {output.uidx2} 2> {log} && ln -s {output.uidx1} {output.idx1} && ln -s {output.uidx2} {output.idx2}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                index2 = rules.generate_index.output.idx2,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -j {input.index2} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                index2 = rules.generate_index.output.idx2,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -j {input.index2} -q {input.query} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"
