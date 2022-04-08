logid = 'LoveSonesonPatro.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
log.info(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config2(SAMPLES, config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
log.info(logid+"COMPARISON: "+str(comparison))

rule themall:
    input:  session = expand("DTU/{combo}/DTU_LSP_SESSION.gz", outdir=outdir),
            results = expand("DTU/{combo}/DTU_LSP_{comparison}_results.tsv.gz", outdir=outdir, comparison=compstr)

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(INDEX),
            uidx = directory(expand("{refd}/INDICES/{mape}_{unikey}", refd=REFDIR, mape=COUNTENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, 'DTU', DTUENV)['OPTIONS'], ['INDEX']))))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  ""+COUNTENV+".yaml"
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
        conda:  ""+COUNTENV+".yaml"
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
        conda:  ""+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'DTU', DTUENV)['OPTIONS'].get('QUANT', ""),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else '-l U',
                linksf = lambda wildcards, output: str(os.path.abspath(output.ctsdir))
        shell: "set +euo pipefail; {params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} &>> {log} && gzip {output.ctsdir}/quant.sf && ln -fs {params.linksf}/quant.sf.gz {output.cnts} &>> {log}"


rule create_annotation_table:
    input:   dir  = expand(rules.mapping.output.ctsdir, file=samplecond(SAMPLES, config)),
    output:  anno = expand("DTU/{combo}/Tables/{scombo}_ANNOTATION.gz", outdir=outdir)
    log:     expand("LOGS/DTU/{combo}/create_DTU_table.log", outdir=outdir)
    conda:   ""+COUNTENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.dir, config,'DTU'),
             bins = BINS
    shell: "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = rules.create_annotation_table.output.anno,
    output: session = rules.themall.input.session,
            results = rules.themall.input.results
    log:    expand("LOGS/DTU/{combo}/run_DTU.log", outdir=outdir)
    conda:  ""+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DTUBIN]),
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            outdir = 'DTU/'+combo,
            ref = ANNOTATION
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.pcombo} {params.compare} {threads} 2> {log}"
