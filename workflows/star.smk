MAPPERBIN, MAPPERENV = env_bin_from_config(config,'MAPPING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = MAPPERENV
unik = get_dict_hash(keydict)

rule generate_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=unik)),
            idxfile = expand("{refd}/INDICES/{mape}_{unikey}/{pref}.idx", refd=REFDIR, mape=MAPPERENV, unikey=unik, pref=PREFIX),
            tmpfa = temp(expand("TMP/{mape}/ref.fa", mape=MAPPERENV)),
            tmpref = temp(expand("TMP/{mape}/ref.anno", mape=MAPPERENV))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda w: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            anno = ANNOTATION,            
            #tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
            pref = PREFIX,
            lnkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmpfa} && zcat {params.anno} > {output.tmpref} && mkdir -p {output.uidx} && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outFileNamePrefix {output.uidx}/{params.pref} --outTmpDir TMP/star_generate_genome --genomeDir {output.uidx} --genomeFastaFiles {output.tmpfa} --sjdbGTFfile {output.tmpref} &> {log} && touch {output.idxfile} && ln -fs {params.lnkidx} {output.idx} && cat {output.uidx}/*Log.out >> {log};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                uidx = rules.generate_index.output.uidx[0],
                dummy = rules.generate_index.output.idxfile,
                ref = REFERENCE
        output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                log = "LOGS/{combo}/{file}/Log.out",
                log_final = "LOGS/{combo}/{file}/Log.final.out",
                tmp = temp("TMP/STAROUT/{combo}/{file}")
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp=MAPPERBIN,
                anno = ANNOTATION,
                pref = PREFIX,
                tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
        shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && gzip -c {output.tmp}.Aligned.out.sam > {output.mapped} && rm -f {output.tmp}.Aligned.out.sam 2>> {log} && touch {output.tmp}.Unmapped.out.mate1 {output.tmp}.Unmapped.out.mate2 && cat {output.tmp}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\" |gzip > {output.unmapped_r1} && cat {output.tmp}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > {output.unmapped_r2} && rm -f {output.tmp}.Unmapped.out* && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp} "

else:
    if paired != 'singlecell':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    uidx = rules.generate_index.output.uidx[0],
                    dummy = rules.generate_index.output.idxfile,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz",
                    log = "LOGS/{combo}/{file}/Log.out",
                    log_final = "LOGS/{combo}/{file}/Log.final.out",
                    tmp = temp("TMP/STAROUT/{combo}/{file}")                    
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN,
                    anno = ANNOTATION,
                    pref = PREFIX,
                    tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
            shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {input.r1} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && gzip -c {output.tmp}.Aligned.out.sam > {output.mapped} && rm -f {output.tmp}.Aligned.out.sam 2>> {log} touch {output.tmp}.Unmapped.out.mate && gzip {output.tmp}.Unmapped.out.mate* && mv {output.tmp}.Unmapped.out.mate*.gz {output.unmapped} 2>> {log} && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp}"
    else:
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    uidx = rules.generate_index.output.uidx[0],
                    dummy = rules.generate_index.output.idxfile,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam.gz", category="MAPPING")),
                    unmapped_r1 = "UNMAPPED/{combo}/{file}_R1_unmapped.fastq.gz",
                    unmapped_r2 = "UNMAPPED/{combo}/{file}_R2_unmapped.fastq.gz",
                    log = "LOGS/{combo}/{file}/Log.out",
                    log_final = "LOGS/{combo}/{file}/Log.final.out",
                    tmp = temp("TMP/STAROUT/{combo}/{file}")
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    stranded = lambda x: '--soloStrand Forward' if stranded == 'fr' else '--soloStrand Reverse' if stranded == 'rf' else '--soloStrand Unstranded',
                    mapp = MAPPERBIN,
                    anno = ANNOTATION,
                    pref = PREFIX,
                    r1 = lambda wildcards, input: input.r1 if not '--soloBarcodeMate 1' in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', "") else input.r2.replace('TRIMMED_','').replace('_trimmed.fastq', '.fastq'),
                    r2 = lambda wildcards, input: input.r2 if not '--soloBarcodeMate 1' in tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', "") else input.r1.replace('TRIMMED_','').replace('_trimmed.fastq', '.fastq').replace({combo}, ''),
                    tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
            shell: "{params.mapp} {params.mpara} {params.stranded} --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN {threads} --genomeDir {input.uidx} --readFilesCommand zcat --readFilesIn {params.r1} {params.r2} --outFileNamePrefix {output.tmp}. --outReadsUnmapped Fastx &> {log} && samtools view -h {output.tmp}.Aligned.sortedByCoord.out.bam | gzip > {output.mapped} && rm -f {output.tmp}.Aligned.sortedByCoord.out.bam; touch {output.tmp}.Unmapped.out.mate1 {output.tmp}.Unmapped.out.mate2 && cat {output.tmp}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\" |gzip > {output.unmapped_r1} && cat {output.tmp}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > {output.unmapped_r2} && rm -f {output.tmp}.Unmapped.out* && mv -f {output.tmp}*Log.out {output.log} && mv -f {output.tmp}*Log.final.out {output.log_final} && rm -f {output.tmp}*.Log.progress.out && mv {output.tmp}*.out* {params.tocopy} 2>> {log} && touch {output.tmp} "
