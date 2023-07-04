MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)

rule generate_index:
    input:  fa = REFERENCE
    output: idx1 = INDEX,
            idx2 = INDEX2,
            uidx1 = expand("{refd}/INDICES/{mape}_{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=unik),
            uidx2 = expand("{refd}/INDICES/{mape}_{unikey}.idx2", refd=REFDIR, mape=MAPPERENV, unikey=unik+'_bs')
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV.replace('bisulfite', '')+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            linkidx1 = lambda wildcards, output: str(os.path.abspath(output.uidx1[0])),
            linkidx2 = lambda wildcards, output: str(os.path.abspath(output.uidx2[0]))
    shell: "{params.indexer} --threads {threads} {params.ipara} -d {input.fa} -x {output.uidx1} -y {output.uidx2} &> {log} && ln -fs {params.linkidx1} {output.idx1} && ln -fs {params.linkidx2} {output.idx2}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx1 = rules.generate_index.output.uidx1[0],
                uidx2 = rules.generate_index.output.uidx2[0],
                fa = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                tmpmap = temp("MAPPED/{combo}/{file}_temp.sam")
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV.replace('bisulfite', '')+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.fa} -i {input.uidx1} -j {input.uidx2} -q {input.r1} -p {input.r2} --threads {threads} -o {output.tmpmap} 2> {log} && cat {output.tmpmap}| tee >(samtools view -h -F 4 |gzip > {output.mapped}) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 {output.unmapped1} -2 {output.unmapped2} - ) 2>> {log} &>/dev/null && touch {output.unmapped1} {output.unmapped2}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                uidx1 = rules.generate_index.output.uidx1[0],
                uidx2 = rules.generate_index.output.uidx2[0],
                fa = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz",
                tmpmap = temp("MAPPED/{combo}/{file}_temp.sam")
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV.replace('bisulfite', '')+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.fa} -i {input.uidx1} -j {input.uidx2} -q {input.query} --threads {threads} -o {output.tmpmap} 2> {log} && cat {output.tmpmap}| tee >(samtools view -h -F 4 |gzip > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 2>> {log} &>/dev/null && touch {output.unmapped}"
