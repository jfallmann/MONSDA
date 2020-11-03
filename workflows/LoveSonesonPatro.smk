DTUBIN, DTUENV = env_bin_from_config3(config,'DTU')
COUNTBIN, COUNTENV = ['salmon','salmon']#env_bin_from_config2(SAMPLES,config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DTU/LoveSonesonPatro"
comparison = comparable_as_string2(config,'DTU')
compstr = [i.split(":")[0] for i in comparison.split(",")]


rule themall:
    input:  expand("COUNTS/Salmon/{file}_counts.sf.gz", file=samplecond(SAMPLES,config)),

rule salmon_index:
    input:  fa = expand("{ref}/{{dir}}/{{trans}}{{name}}.fa.gz", ref=REFERENCE),
    output: idx = expand("{ref}/{{dir}}/{map}/{{extension}}/{{trans}}{{name}}_{{extension}}/{map}.idx", ref=REFERENCE, map=COUNTENV),
    log:    expand("LOGS/{{dir}}/{{trans}}{{name}}_{{extension}}_{map}.idx.log", map=COUNTENV)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'DTU')['OPTIONS'][0].items()),
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(transcriptomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'DTU')['ANNOTATION']]),
            transpath = lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/{trans}{name}_{extension}".format(ref=REFERENCE, dir=wildcards.dir, trans=wildcards.trans, name=wildcards.name, map=COUNTENV, extension=check_tool_params(SAMPLES[0], None ,config, 'DTU',2))),
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP'),
    shell:  "rm -rf {params.tmpidx} && if [[ -f \"{params.transpath}\" ]]; then ln -fs {params.transpath} {output.idx} && echo \"Found Salmon index, continue with quantify\" ; else {params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.idx} 2>> {log};fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.salmon_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), trans=transcriptome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=COUNTENV, extension=check_tool_params(wildcards.file, None ,config, 'DTU',2))
        output: cnts = report("COUNTS/Salmon/{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/Salmon/{file}", category="COUNTING")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.salmon_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), trans=transcriptome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=COUNTENV, extension=check_tool_params(wildcards.file, None ,config, 'DTU',2))
        output: cnts = report("COUNTS/Salmon/{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/Salmon/{file}", category="COUNTING")
        log:    "LOGS/{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"

# rule prepare_count_table:
#     input:   cnd  = expand(rules.mapping.output.cnts, file=samplecond(SAMPLES,config))
#     output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
#              anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
#     log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
#     conda:   "nextsnakes/envs/"+DTUENV+".yaml"
#     threads: 1
#     params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DTU'),
#              bins = BINS,
#              tpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(samplecond(SAMPLES,config)[0], None ,config, "DTU")['OPTIONS'][1].items())
#     shell: "{params.bins}/Analysis/build_salmon_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} {params.tpara} --loglevel DEBUG 2> {log}"
#
# rule run_dexseqDTU:
#     input:  tbl = rules.prepare_count_table.output.tbl,
#             anno = rules.prepare_count_table.output.anno,
#     output: all = rules.themall.input.all,
#             sum = rules.themall.input.allsum,
#             tbl = rules.themall.input.tbl,
#             bcv = rules.themall.input.bcv,
#             qld = rules.themall.input.qld,
#             dift = rules.themall.input.dift,
#             plot = rules.themall.input.plot,
#             session = rules.themall.input.session
#     log:    expand("LOGS/{outdir}run_dexseq.log",outdir=outdir)
#     conda:  "nextsnakes/envs/"+DTUENV+".yaml"
#     threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
#     params: bins   = str.join(os.sep,[BINS,DTUBIN]),
#             outdir = outdir,
#             compare = comparison
#     shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.outdir} {params.compare} {threads} 2> {log}"
