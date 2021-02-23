MAPPERBIN, MAPPERENV = env_bin_from_config3(config, 'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'][0]))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'][0].items()),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell: "{params.indexer} --threads {threads} {params.ipara} -d {input.fa} -x {output.uidx} 2> {log} && ln -fs {params.linkidx} {output.idx}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -q {input.query} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"
