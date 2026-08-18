COUNTBIN, COUNTENV = env_bin_from_config(config,'COUNTING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = COUNTENV
unik = get_dict_hash(keydict)

bammode = tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('BAM', "").lower()

rule themall:
    input:  expand("COUNTS/{combo}/{file}_counts.sf.gz", combo=combo, file=samplecond(SAMPLES, config))

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=unik))
    log:    expand("LOGS/{sets}/COUNTING/{cape}/idx.log", sets=SETS, cape=COUNTENV)
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('INDEX', ""),
            decoy = f"-d {os.path.abspath(DECOY)}" if DECOY else '',
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} {params.decoy} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"


if bammode == 'transcriptome':
    rule mapping:
        input:  bam = "MAPPED/{combo}/{file}_mapped_sorted.bam",
                fa = REFERENCE
        output: cnts = report("COUNTS/{combo}/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/COUNTING/salmon/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                nbam = lambda wildcards, output: os.path.join(str(output.ctsdir), "namecollated.bam"),
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; mkdir -p {output.ctsdir}; samtools collate -@ {threads} -o {params.nbam} {input.bam} &>> {log}; {params.mapp} quant -p {threads} -t {input.fa} -a {params.nbam} --deterministic {params.cpara} -o {output.ctsdir} &>> {log} && rm -f {params.nbam} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

elif bammode == 'genome':
    rule mapping:
        input:  bam = "MAPPED/{combo}/{file}_mapped_sorted.bam",
                anno = ANNOTATION,
                fa = REFERENCE
        output: cnts = report("COUNTS/{combo}/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/COUNTING/salmon/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                nbam = lambda wildcards, output: os.path.join(str(output.ctsdir), "namecollated.bam"),
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; mkdir -p {output.ctsdir}; samtools collate -@ {threads} -o {params.nbam} {input.bam} &>> {log}; {params.mapp} quant -p {threads} -a {params.nbam} --annotation {input.anno} --genome {input.fa} {params.cpara} -o {output.ctsdir} &>> {log} && rm -f {params.nbam} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

elif paired == 'paired':
                r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo),
                uidx = rules.salmon_index.output.uidx[0]
        output: cnts = report("COUNTS/{combo}/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/COUNTING/salmon/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else '-l IU',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.uidx} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

else:
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo),
                uidx = rules.salmon_index.output.uidx[0]
        output: cnts = report("COUNTS/{combo}/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/COUNTING/salmon/salmonquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.uidx} {params.stranded} {params.cpara} -o {output.ctsdir} -r {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf ; ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"
