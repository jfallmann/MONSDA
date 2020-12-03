logid = 'drimseq_DTU.smk '
DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
log.info(logid+"DTUENV: "+str(DTUENV))
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config2(SAMPLES,config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DTU/DRIMSEQ"
comparison = comparable_as_string2(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]
log.info(logid+"COMPARISON: "+str(comparison))

rule themall:
    input:  session = expand("{outdir}/DRIMSEQ_DTU_SESSION.gz", outdir=outdir),
            res1 = expand("{outdir}/DTU_DRIMSEQ_{comparison}_results.tsv.gz", outdir=outdir, comparison=compstr),
            res2 = expand("{outdir}/DTU_DRIMSEQ_{comparison}_stageR-filtered.tsv.gz", outdir=outdir, comparison=compstr),
            res3 = expand("{outdir}/DTU_DRIMSEQ_{comparison}_post-hoc-filtered-on-SD.tsv.gz", outdir=outdir, comparison=compstr)

rule salmon_index:
    input:  fa = REFERENCE
    output: idx = directory(expand("{refd}/INDICES/{mape}", refd=REFDIR, mape=COUNTENV))
    log:    expand("LOGS/{sets}/salmon.idx.log", sets=SETS)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'DTU')['OPTIONS'][0].items()),
    shell:  "{params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.idx} 2>> {log}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "FASTQ/{file}_R1.fastq.gz",
                r2 = "FASTQ/{file}_R2.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "FASTQ/{file}.fastq.gz",
                index = rules.salmon_index.output.idx
        output: ctsdir = directory("COUNTS/Salmon/{file}")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} "

rule create_annotation_table:
    input:   dir  = expand(rules.mapping.output.ctsdir, file=samplecond(SAMPLES,config)),
    output:  anno = expand("{outdir}/Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}/create_DTU_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.dir,config,'DTU'),
             bins = BINS
             # tpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(samplecond(SAMPLES,config)[0], None ,config, "DTU")['OPTIONS'][1].items())
    shell: "python3 {params.bins}/Analysis/build_DTU_table.py {params.dereps} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_DTU:
    input:  anno = rules.create_annotation_table.output.anno,
    output: session = rules.themall.input.session,
            res1 = rules.themall.input.res1,
            res2 = rules.themall.input.res2,
            res3 = rules.themall.input.res3
    log:    expand("LOGS/{outdir}run_DTU.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DTUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DTUBIN]),
            compare = comparison,
            outdir = outdir,
            ref = ANNOTATION,
            # pvcut = lambda wildcards: ' '.join(f"{val}" for (key,val) in tool_params(SAMPLES[0], None ,config, 'DTU')['OPTIONS'][2].items())
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {params.ref} {params.outdir} {params.compare} {threads} 2> {log}"
