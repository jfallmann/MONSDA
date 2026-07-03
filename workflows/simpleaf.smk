COUNTBIN, COUNTENV = env_bin_from_config(config,'COUNTING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = COUNTENV
unik = get_dict_hash(keydict)

rule themall:
    input:  expand("COUNTS/{combo}/{file}/af_quant/alevin/quants_mat.mtx.gz", combo=combo, file=samplecond(SAMPLES, config))

rule generate_index:
    input:  fa = REFERENCE,
            anno = ANNOTATION
    output: idx = directory(expand("{refd}/INDICES/{cenv}_{unikey}", refd=REFDIR, cenv=COUNTENV, unikey=unik)),
            idxfile = expand("{refd}/INDICES/{cenv}_{unikey}/index/info.json", refd=REFDIR, cenv=COUNTENV, unikey=unik),
            tmpfa = temp(expand("TMP/{cenv}/ref.fa", cenv=COUNTENV)),
            tmpanno = temp(expand("TMP/{cenv}/ref.gtf", cenv=COUNTENV))
    log:    expand("LOGS/{sets}/{cenv}.idx.log", sets=SETS, cenv=COUNTENV)
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: mapp = COUNTBIN,
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('INDEX', "--ref-type spliced+unspliced"),
            afhome = lambda wildcards, output: os.path.join(os.path.dirname(str(output.idx[0])), "af_home")
    shell:  "if [[ -f \"{output.idxfile}\" ]]; then touch {output.idxfile} && echo \"Found index, continue with quant\" ; else export ALEVIN_FRY_HOME={params.afhome}; mkdir -p {params.afhome}; {params.mapp} set-paths &>> {log}; zcat {input.fa} > {output.tmpfa} && zcat {input.anno} > {output.tmpanno} && {params.mapp} index --fasta {output.tmpfa} --gtf {output.tmpanno} -t {threads} {params.ipara} -o {output.idx} &> {log};fi"

if paired == 'paired' or paired == 'singlecell':
    rule counting:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                idx = rules.generate_index.output.idx,
                dummy = rules.generate_index.output.idxfile
        output: cnts = report("COUNTS/{combo}/{file}/af_quant/alevin/quants_mat.mtx.gz", category="COUNTING"),
                ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
        log:    "LOGS/{combo}/{file}/simpleafquant.log"
        conda:  ""+COUNTENV+".yaml"
        container: "oras://jfallmann/monsda:"+COUNTENV+""
        threads: MAXTHREAD
        params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', "--chemistry 10xv3 --resolution cr-like --expected-ori fw --unfiltered-pl"),
                mapp = COUNTBIN,
                afhome = lambda wildcards, input: os.path.join(os.path.dirname(str(input.idx)), "af_home"),
                quantdir = lambda wildcards, output: os.path.join(str(output.ctsdir), "af_quant")
        shell: "export ALEVIN_FRY_HOME={params.afhome}; mkdir -p {output.ctsdir}; {params.mapp} quant --index {input.idx}/index --reads1 {input.r1} --reads2 {input.r2} -t {threads} {params.cpara} -o {params.quantdir} &> {log} && gzip -c {params.quantdir}/alevin/quants_mat.mtx > {output.cnts} 2>> {log}"
