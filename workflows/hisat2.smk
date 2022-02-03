MAPPERBIN, MAPPERENV = env_bin_from_config3(config,'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])))),
            idxfile = expand("{refd}/INDICES/{mape}_{unikey}/{pref}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])), pref = PREFIX),
            tmp = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0])),
            pref = PREFIX
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then ln -fs {params.lnkidx} {output.idx} && mkdir -p {output.uidx} && touch {output.tmp} {output.idxfile} && echo \"Found hisat index, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && mkdir -p {output.uidx} ;{params.indexer}-build {params.ipara} -p {threads} {output.tmp} {output.uidx}/{params.pref} 2> {log} && ln -fs {params.lnkidx} {output.idx} && touch {output.uidx} && touch {output.idxfile};fi"

hs2alg = MAPPERBIN.split(' ')[1] if ' ' in MAPPERBIN else MAPPERBIN

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                summary = "MAPPED/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else '',
                um = lambda wildcards, output: output.unmapped_r1.replace('_R1_unmapped.fastq.gz', '.unmapped.gz'),
                umn = lambda wildcards, output: output.unmapped_r1.replace('_R1_unmapped.fastq.gz', ''),
                pref = PREFIX
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {input.uidx}/{params.pref} -1 {input.r1} -2 {input.r2} -S {output.mapped} --un-conc-gz {params.um} --new-summary --summary-file {output.summary} 2>> {log} && touch {params.umn}.unmapped.1.gz {params.umn}.unmapped.2.gz; rename 's/.unmapped.([1|2]).gz/_R$1_unmapped.fastq.gz/' {params.umn}.unmapped.*.gz; touch {output.unmapped_r1} {output.unmapped_r2} &>> {log}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz",
                summary = "MAPPED/{combo}/{file}.summary"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else '',
                pref = PREFIX
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {input.uidx}/{params.pref} -U {input.query} -S {output.mapped} --un-gz {output.unmapped} --new-summary --summary-file {output.summary} 2>> {log} && touch {output.unmapped}"
