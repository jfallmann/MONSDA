DTUBIN, DTUENV = env_bin_from_config(config,'DTU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DTU/DEXSEQ"
comparison = comparable_as_string(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
keydict = subDict(tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
unik = get_dict_hash(keydict)

rule themall:
    input:


rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=unik))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  "MONSDA/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'].get('INDEX', ""),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0]))
    shell:  "set +euo pipefail; {params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.uidx} &>> {log} && ln -fs {params.linkidx} {output.idx}"


if paired == 'paired':
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R1_trimmed.fastq.gz", scombo=scombo) if not rundedup else expand("DEDUP_FASTQ/{scombo}/{{file}}_R1_dedup.fastq.gz", scombo=scombo),
                r2 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_R2_trimmed.fastq.gz", scombo=scombo) if not rundedup else expand("DEDUP_FASTQ/{scombo}/{{file}}_R2_dedup.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "MONSDA/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else '-l IU',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"

else:
    rule mapping:
        input:  r1 = expand("TRIMMED_FASTQ/{scombo}/{{file}}_trimmed.fastq.gz", scombo=scombo) if not rundedup else expand("DEDUP_FASTQ/{scombo}/{{file}}_dedup.fastq.gz", scombo=scombo),
                index = rules.salmon_index.output.idx,
                uix = rules.salmon_index.output.uidx
        output: cnts = report("DTU/{combo}/salmon/{file}_counts.sf.gz", category="COUNTING"),
                ctsdir = report(directory("DTU/{combo}/salmon/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/salmonquant.log"
        conda:  "MONSDA/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"
