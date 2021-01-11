COUNTBIN, COUNTENV = env_bin_from_config3(config,'COUNTING')

rule themall:
    input:  expand("COUNTS/Salmon/{combo}{file}_counts.sf.gz", file=samplecond(SAMPLES,config)),

rule salmon_index:
    input:  fa = REFERENCE,
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}/{unikey}.idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))
    log:    expand("LOGS/{sets}/{cape}.idx.log", sets=SETS, cape=COUNTENV)
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'][0].items()),
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(transcriptomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['ANNOTATION']]),
            transpath = expand("{refd}/{mape}/{unikey}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'][0])),  #lambda wildcards: os.path.abspath("{ref}/{dir}/{map}/{extension}/{trans}{name}_{extension}".format(ref=REFERENCE, dir=wildcards.dir, trans=wildcards.trans, name=wildcards.name, map=COUNTENV, extension=check_tool_params(SAMPLES[0], None ,config, 'COUNTING',2))),
            linkidx = lambda wildcards, output: str(os.path.abspath(output.uidx[0])),
            tmpidx = lambda x: tempfile.mkdtemp(dir='TMP')
    shell:  "rm -rf {params.tmpidx} && if [[ -f \"{params.transpath}\" ]]; then ln -fs {params.transpath} {output.idx} && echo \"Found Salmon index, continue with quantify\" ; else {params.mapp} index {params.ipara} -p {threads} -t {input.fa} -i {output.uidx} 2>> {log} && ln -fs {params.linkidx} {output.idx} ;fi"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}{file}_R2_trimmed.fastq.gz",
                index = rules.generate_index.output.idx
        output: cnts = report("COUNTS/Salmon/{combo}{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/Salmon/{combo}{file}", category="COUNTING")
        log:    "LOGS/{combo}{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'COUNTING', COUNTENV)['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l ISF' if (stranded == 'fr' or stranded == 'ISF') else '-l ISR' if (stranded == 'rf' or stranded == 'ISR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} -2 {input.r2} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"

else:
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}{file}_trimmed.fastq.gz",
                index = rules.generate_index.output.idx
        output: cnts = report("COUNTS/Salmon/{combo}{file}.sf.gz", category="COUNTING"),
                ctsdir = report("COUNTS/Salmon/{combo}{file}", category="COUNTING")
        log:    "LOGS/{combo}{file}/salmonquant.log"
        conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
        threads: MAXTHREAD
        params: cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'COUNTING', COUNTENV)['OPTIONS'][1].items()),
                mapp=COUNTBIN,
                stranded = lambda x: '-l SF' if (stranded == 'fr' or stranded == 'SF') else '-l SR' if (stranded == 'rf' or stranded == 'SR') else ''
        shell: "{params.mapp} quant -p {threads} -i {input.index} {params.stranded} {params.cpara} -o {output.ctsdir} -1 {input.r1} 2>> {log} && gzip {output.ctsdir}/quant.sf && ln -s {output.ctsdir}/quant.sf.gz {output.cnts} 2>> {log}"
