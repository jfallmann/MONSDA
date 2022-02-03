MAPPERBIN, MAPPERENV = env_bin_from_config3(config,'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])))),
            dummy = expand("{refd}/INDICES/{mape}_{unikey}/{pref}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])), pref=PREFIX)
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda w: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            anno = ANNOTATION,            
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
            pref = PREFIX,
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "if [[ -f \"{output.dummy}\" ]]; then touch {output.dummy} && ln -fs {output.dummy} {output.idx} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {params.tmpidx}/star_ref.fa && zcat {params.anno} > {params.tmpidx}/star_ref.anno && mkdir -p {output.uidx} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {output.uidx}/{params.pref} --outTmpDir {params.tmpidx}/star --genomeDir {output.uidx} --genomeFastaFiles {params.tmpidx}/star_ref.fa --sjdbGTFfile {params.tmpidx}/star_ref.anno 2> {log} && touch {output.dummy} && ln -fs {params.lnkidx} {output.idx} && cat {output.uidx}/{params.pref}Log.out >> {log} && rm -rf {params.tmpidx};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.dummy,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{combo}/{file}_unmapped_R1.fastq.gz",
                unmapped_r2 = "UNMAPPED/{combo}/{file}_unmapped_R2.fastq.gz",
                tmp = temp("TMP/STAROUT/{combo}/{file}")
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                anno = ANNOTATION,
                pref = PREFIX,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate1 && mv {output.tmp}.Unmapped.out.mate1.gz {output.unmapped_r1} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate2 && mv {output.tmp}.Unmapped.out.mate2.gz {output.unmapped_r2} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.dummy,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz",
                tmp = temp("TMP/STAROUT/{combo}/{file}")
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                anno = ANNOTATION,
                pref = PREFIX,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate* && mv {output.tmp}.Unmapped.out.mate*.gz {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
