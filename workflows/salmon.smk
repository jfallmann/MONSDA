COUNTBIN, COUNTENV = env_bin_from_config3(config,'COUNTING')

rule themall:
    input:  expand("COUNTS/{combo}/{file}_counts.sf.gz", file=samplecond(SAMPLES, config)),

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=COUNTENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'][0]))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  "NextSnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'][0].items()),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "{params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.uidx} 2>> {log} && ln -fs {params.linkidx} {output.idx}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{scombo}/{file}_R1_trimmed.fastq.gz" if not rundedup else "DEDUP_FASTQ/{scombo}/{file}_R1_dedup.fastq.gz",
                r2 = "TRIMMED_FASTQ/{scombo}/{file}_R2_trimmed.fastq.gz" if not rundedup else "DEDUP_FASTQ/{scombo}/{file}_R2_dedup.fastq.gz",
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("COUNTS/{combo}/{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/Salmon/{combo}/{file}", category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "NextSnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{scombo}/{file}_trimmed.fastq.gz" if not rundedup else "DEDUP_FASTQ/{scombo}/{file}_dedup.fastq.gz",
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("COUNTS/{combo}/{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/{combo}/{file}", category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "NextSnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"
