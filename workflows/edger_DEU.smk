DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DEU/EDGER/"
comparison=comparable_as_string(config,'DEU')

rule themall:
    input:  all = expand("{outdir}All_Conditions_MDS.png", outdir=outdir),
            tbl = expand("{outdir}normalized_table_{comparison}.tsv", outdir=outdir, comparison=comparison.split(",")),
            plot = expand("{outdir}{comparison}_MD.png", outdir=outdir, comparison=comparison.split(",")),
            bcv = expand("{outdir}{comparison}_BCV.png", outdir=outdir, comparison=comparison.split(",")),
            qld = expand("{outdir}{comparison}_QLDisp.png", outdir=outdir, comparison=comparison.split(",")),
            mds = expand("{outdir}{comparison}_MDS.png", outdir=outdir, comparison=comparison.split(","))

rule featurecount_unique:
    input:  mapf = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts = "COUNTS/Featurecounter_edger_deu/{file}_mapped_sorted_unique.counts",
            anno = "COUNTS/Featurecounter_edger_deu/{file}_edger.gtf.gz"
    log:    "LOGS/{file}/deu_edger_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count  = COUNTBIN,
            anno   = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]),
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            countstrand = lambda x: '-s' if stranded == 'fr' or stranded == 'rf' else '',
            countgtf = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_fc_edger.gtf')),
            dexgtf   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_edger.gtf'))
    shell:  "if [ ! -f \"{params.dexgtf}\" ] || [ ! -f \"{params.countgtf}\" ];then {params.bins}/Analysis/DEU/prepare_dexseq_annotation2.py -f {params.countgtf} {params.countstrand} {params.anno} {params.dexgtf} ;fi && ln -s {params.dexgtf} {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.countgtf}) -o {output.cts} {input.mapf} 2> {log}"


rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = "DEU/Tables/EDGER/RUN_DEU_Analysis.tbl.gz",
             anno = "DEU/Tables/EDGER/RUN_DEU_Analysis.anno.gz"
    log:     "LOGS/DEU/prepare_count_table.log"
    conda:   "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS
    shell: "{params.bins}/Analysis/DEU/build_DEXSeq_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_edger:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.themall.input.all,
            rules.themall.input.tbl,
            rules.themall.input.plot,
            rules.themall.input.bcv,
            rules.themall.input.qld,
            rules.themall.input.mds
    log:    "LOGS/DEU/run_edger.log"
    conda:  "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison,
            flat   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_edger.gtf'))
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.flat} {params.outdir} {params.compare} {threads} 2> {log} "


onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
