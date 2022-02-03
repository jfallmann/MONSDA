MAPPERBIN, MAPPERENV = env_bin_from_config3(config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx1 = INDEX,
	        idx2 = INDEX2,
            uidx1 = expand("{refd}/INDICES/{mape}_{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subdict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])))
            uidx2 = expand("{refd}/INDICES/{mape}_{unikey}.idx2", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])))+'second'
    log:    expand("LOGS/{mape}.idx.log", mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'].get('INDEX', ""),
            linkidx1 = lambda wildcards, output: str(os.path.abspath(output.uidx1[0])),
            linkidx2 = lambda wildcards, output: str(os.path.abspath(output.uidx2[0]))
    shell: "{params.indexer} --threads {threads} {params.ipara} -d {input.fa} -x {output.uidx1} -y {output.uidx2} &> {log} && ln -fs {params.linkidx1} {output.idx1} && ln -s {params.linkidx2} {output.idx2}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx1 = rules.generate_index.output.uidx[0],
                uidx2= rules.generate_index.output.uidx2[0],
                ref = REFERENCE
        output: mapped = report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.uidx1} -j {input.uidx2} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} &> {log}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                uidx1 = rules.generate_index.output.uidx[0],
                uidx2= rules.generate_index.output.uidx2[0],
                ref = REFERENCE
        output: mapped = report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.uidx1} -j {input.uidx2} -q {input.query} --threads {threads} -o {output.mapped} -u {output.unmapped} &> {log}"
