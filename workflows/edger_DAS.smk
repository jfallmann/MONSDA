DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/EDGER/"
comparison=comparable_as_string(config,'DAS')

rule themall:
    input:  all = expand("{outdir}All_Conditions_MDS.png", outdir=outdir),
            tbl = expand("{outdir}{comparison}_normalized_table_.tsv", outdir=outdir, comparison=comparison.split(",")),
            plot = expand("{outdir}{comparison}_MD.png", outdir=outdir, comparison=comparison.split(",")),
            bcv = expand("{outdir}{comparison}_BCV.png", outdir=outdir, comparison=comparison.split(",")),
            qld = expand("{outdir}{comparison}_QLDisp.png", outdir=outdir, comparison=comparison.split(",")),
            dift = expand("{outdir}{comparison}_diffSplice_{test}.png", outdir=outdir, comparison=comparison.split(","), test=["geneTest","simesTest","exonTest"])
            tops = expand("{outdir}{comparison}_topSplice_simes_{n}.png", outdir=outdir, comparison=comparison.split(","), n=range(1,11))

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = "COUNTS/Featurecounter_edger/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecount_edger_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DAS")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = "DAS/Tables/EDGER/COUNTS.gz",
             anno = "DAS/Tables/EDGER/ANNOTATION.gz"
    log:     "LOGS/DAS/prepare_count_table.log"
    conda:   "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DAS'),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edger:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.themall.input.all,
            rules.themall.input.tbl,
            rules.themall.input.plot,
            rules.themall.input.bcv,
            rules.themall.input.qld,
            rules.themall.input.dift,
            rules.themall.input.tops
    log:    "LOGS/DAS/run_edger.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.outdir} {params.compare} {threads} 2> {log} "

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
