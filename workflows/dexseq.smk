DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DEU/DEXSEQ/"
comparison=comparable_as_string(config,'DEU')

rule themall:
    input: tbl = expand("{outdir}DEXSeq_{comparison}.tsv.gz", outdir=outdir, comparison=comparison.split(",")),
           plot = expand("{outdir}DEXSeq_{comparison}_DispEsts.pdf", outdir=outdir, comparison=comparison.split(",")),
           html = expand("{outdir}DEXSeqReport_{comparison}/DEXSeq_{comparison}.html", outdir=outdir, comparison=comparison.split(","))
           session = expand("{outdir}DEXSeq_SESSION.gz", outdir=outdir)# R object?

rule featurecount_dexseq_unique:
    input:  mapf = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts  = "COUNTS/Featurecounter_dexseq/{file}_mapped_sorted_unique.counts",
            anno = "COUNTS/Featurecounter_dexseq/{file}_dexseq.gtf.gz"
    log:    "LOGS/{file}/featurecount_dexseq_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count  = COUNTBIN,
            anno   = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]),
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            countstrand = lambda x: '-s' if stranded == 'fr' or stranded == 'rf' else '',
            countgtf = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_fc_dexseq.gtf')),
            dexgtf   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_dexseq.gtf'))
    shell:  "if [ ! -f \"{params.dexgtf}\" ] || [ ! -f \"{params.countgtf}\" ];then {params.bins}/Analysis/DEU/prepare_dexseq_annotation2.py -f {params.countgtf} {params.countstrand} {params.anno} {params.dexgtf} ;fi && ln -s {params.dexgtf} {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.countgtf}) -o {output.cts} {input.mapf} 2> {log}"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_dexseq_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = "DEU/Tables/RUN_DEU_Analysis.tbl.gz",
             anno = "DEU/Tables/RUN_DEU_Analysis.anno.gz"
    log:     "LOGS/DEU/prepare_count_table.log"
    conda:   "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS
    shell: "{params.bins}/Analysis/DEU/build_DEXSeq_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_dexseq:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno
    output: plot = rules.themall.input.plot,
            tbl  = rules.themall.input.tbl,
            html = rules.themall.input.html,
            session = rules.themall.input.session
    log:    "LOGS/DEU/run_dexseq.log"
    conda:  "snakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison,
            flat   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_dexseq.gtf')),
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {params.flat} {params.outdir} {params.compare} {threads} 2> {log}"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
