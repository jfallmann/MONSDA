MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)

rule generate_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}_{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=unik)
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer=MAPPERBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell: "{params.indexer} -t {threads} -d {output.uidx} {params.ipara} {input.fa} 2> {log} && ln -s {params.linkidx} {output.idx}"

if paired == 'paired':
    rule mapping:
        input:  q1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                q2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0]
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN
        shell: "{params.mapp} -t {threads} {params.mpara} {input.uidx} {input.q1} {input.q2} 2> {log}| tee >(samtools view -h -F 4 |gzip > {output.mapped}) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 {output.unmapped1} -2 {output.unmapped2} - ) 2>> {log} &>/dev/null && touch {output.unmapped1} {output.unmapped2}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0]
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN
        shell: "{params.mapp} -t {threads} {params.mpara} {input.uidx} {input.query} 2> {log}| tee >(samtools view -h -F 4 |gzip > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 2>> {log} &>/dev/null && touch {output.unmapped}"
