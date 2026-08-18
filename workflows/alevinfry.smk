COUNTBIN, COUNTENV = env_bin_from_config(config,'COUNTING')
keydict = sub_dict(tool_params(SAMPLES[0], None, config, 'COUNTING', COUNTENV)['OPTIONS'], ['INDEX'])
keydict["REF"] = REFERENCE
keydict["DECOY"] = DECOY
keydict["ENV"] = COUNTENV
unik = get_dict_hash(keydict)

rule themall:
    input:  expand("COUNTS/{combo}/{file}/alevin/quants_mat.mtx.gz", combo=combo, file=samplecond(SAMPLES, config))

rule generate_t2g:
    input:  anno = ANNOTATION
    output: t2g = expand("{refd}/INDICES/{cenv}_{unikey}/t2g.tsv", refd=REFDIR, cenv=COUNTENV, unikey=unik)
    log:    expand("LOGS/{sets}/COUNTING/{cenv}/t2g.log", sets=SETS, cenv=COUNTENV)
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: 1
    shell:  "mkdir -p $(dirname {output.t2g}); zcat {input.anno} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}} $3==\"transcript\"{{match($9,/transcript_id \"[^\"]+\"/); t=substr($9,RSTART+15,RLENGTH-16); match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10); print t,g}}' > {output.t2g} 2> {log}"

rule counting:
    input:  rad = "MAPPED/{combo}/{file}",
            radfile = "MAPPED/{combo}/{file}/map.rad",
            t2g = rules.generate_t2g.output.t2g
    output: cnts = report("COUNTS/{combo}/{file}/alevin/quants_mat.mtx.gz", category="COUNTING"),
            ctsdir = report(directory("COUNTS/{combo}/{file}"), category="COUNTING")
    log:    "LOGS/{combo}/{file}/COUNTING/alevinfry/alevinfryquant.log"
    conda:  ""+COUNTENV+".yaml"
    container: "oras://jfallmann/monsda:"+COUNTENV+""
    threads: MAXTHREAD
    params: cpara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('COUNT', ""),
            ppara = lambda wildcards: tool_params(wildcards.file, None, config, 'COUNTING', COUNTENV)['OPTIONS'].get('PERMIT', "--unfiltered-pl"),
            mapp = COUNTBIN,
            quantdir = lambda wildcards, output: os.path.join(str(output.ctsdir), "af_quant")
    shell: "{params.mapp} generate-permit-list -i {input.rad} -d fw {params.ppara} -o {params.quantdir} &> {log}; {params.mapp} collate -i {params.quantdir} -r {input.rad} -t {threads} &>> {log}; {params.mapp} quant -i {params.quantdir} -m {input.t2g} -t {threads} -r cr-like {params.cpara} -o {params.quantdir} --use-mtx &>> {log} && mkdir -p $(dirname {output.cnts}) && gzip -c {params.quantdir}/alevin/quants_mat.mtx > {output.cnts} 2>> {log}"
