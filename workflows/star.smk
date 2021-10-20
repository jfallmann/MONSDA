MAPPERBIN, MAPPERENV = env_bin_from_config3(config,'MAPPING')

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = expand("{refd}/INDICES/{mape}_{unikey}/{pref}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])), pref=PREFIX)
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda w: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            anno = ANNOTATION,
            linkidx = lambda wildcards, output: str(os.path.abspath(str.join(os.sep, str(output.uidx[0]).split(os.sep)[:-1]))) if PREFIX != '' else str(os.path.abspath(str(output.uidx[0]))),
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
            pref = PREFIX
    shell:  "if [[ -f \"{params.linkidx}/{params.pref}SAindex\" ]]; then ln -fs {params.linkidx} {output.idx} && touch {output.uidx} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {params.tmpidx}/star_ref.fa && zcat {params.anno} > {params.tmpidx}/star_ref.anno && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {params.linkidx}/{params.pref} --outTmpDir {params.tmpidx}/star --genomeDir {params.linkidx} --genomeFastaFiles {params.tmpidx}/star_ref.fa --sjdbGTFfile {params.tmpidx}/star_ref.anno 2> {log} && ln -fs {params.linkidx} {output.idx} && touch {output.uidx} && cat {params.linkidx}/{params.pref}Log.out >> {log} && rm -rf {params.tmpidx};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                uidx = rules.generate_index.output.uidx,
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
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.index} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate1 && mv {output.tmp}.Unmapped.out.mate1.gz {output.unmapped_r1} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate2 && mv {output.tmp}.Unmapped.out.mate2.gz {output.unmapped_r2} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx,
                uidx = rules.generate_index.output.uidx,
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
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.index} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx 2>> {log} && mv {output.tmp}.Aligned.out.sam {output.mapped} 2>> {log} && gzip {output.tmp}.Unmapped.out.mate* && mv {output.tmp}.Unmapped.out.mate*.gz {output.unmapped} 2>> {log} && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
