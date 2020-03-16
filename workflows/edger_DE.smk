DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DE/EDGER/"
comparison=comparable_as_string(config,'DE')

rule themall:
    input:  plot = expand("{outdir}{comparison}.png", outdir=outdir, comparison=comparison.split(",")),
            bcv = expand("{outdir}BCV.png", outdir=outdir),
            mds = expand("{outdir}MDS.png", outdir=outdir),
            tbl = expand("{outdir}normalized_table.tsv", outdir=outdir)

rule featurecount_unique:
    input:  mapf = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts = "COUNTS/Featurecounter_edger_de/{file}_mapped_sorted_unique.counts"
            #temp("COUNTS/Featurecounter_edger_de/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_de_edger_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DE')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DE")['OPTIONS'][0].items()),
            #anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            #cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t gene -g '+config['COUNTING']['FEATURES']['gene'],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output[0]} {input[0]} 2> {log}"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/EDGER/RUN_DE_Analysis.tbl.gz",
             anno = "DE/Tables/EDGER/RUN_DE_Analysis.anno.gz"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.cnd,config,'DE'),
            #decond = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')["GROUP"]) for x in input.cnd]),
            #samples = lambda wildcards, input: str.join(',',input.cnd),
            bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_edger:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.themall.input.plot,
            rules.themall.input.bcv,
            rules.themall.input.mds,
            rules.themall.input.tbl
    log:    "LOGS/DE/run_edger.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params: bins   = str.join(os.sep,[BINS,DEBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.tbl} {input.anno} {params.outdir} {params.compare} 2> {log} "


onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
