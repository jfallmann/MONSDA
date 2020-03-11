DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DE/DESEQ2/"
comparison=comparable_as_string(config,'DE')

rule all:
    input:  plot = expand("{outdir}{comparison}_DESeq2_plot.pdf", outdir=outdir, comparison=comparison.split(",")),
            rld = expand("{outdir}{comparison}_DESeq2_rld.txt.gz", outdir=outdir, comparison=comparison.split(",")),
            vsd = expand("{outdir}{comparison}_DESeq2_vsd.txt.gz", outdir=outdir, comparison=comparison.split(",")),
            csv = expand("{outdir}{comparison}.csv.gz", outdir=outdir, comparison=comparison.split(",")),
            heat = expand("{outdir}DESeq2_heatmap{i}.pdf", outdir=outdir,i=[1,2,3,"_samplebysample"]),
            pca = expand("{outdir}DESeq2_PCA.pdf", outdir=outdir),
            vst = expand("{outdir}DESeq2_VST_and_log2.pdf", outdir=outdir),
            rpl = expand("{outdir}Rplots.pdf", outdir=outdir),

rule featurecount_unique:
    input:  reads = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: cts   = "COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts",
            anno  = temp("COUNTS/Featurecounter_genes/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_de_gene_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t gene -g '+config['COUNTING']['FEATURES']['gene'],
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output.anno} -o {output.cts} {input.reads} 2> {log}"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = "DE/Tables/RUN_DE_Analysis.tbl.gz",
             anno = "DE/Tables/RUN_DE_Analysis.anno.gz"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  decond = lambda wildcards, input: str.join(',',get_reps(input.cnd,config,'DE','CONDITIONS')),
             dereps = lambda wildcards, input: str.join(',',get_reps(input.cnd,config,'DE','REPLICATES')),
             detypes = lambda wildcards, input: '-t '+str.join(',',get_reps(input.cnd,config,'DE','TYPES')),
             paired = lambda wildcards, input:  '--paired '+str.join(',',[checkpaired_rep([str.join(os.sep,x.split(os.sep)[2:]) for x in get_reps(input.cnd,config,'DE','REPLICATES')],config)]),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py -r {params.dereps} -c {params.decond} {params.detypes} {params.paired} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_deseq2:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.all.input.plot,
            rules.all.input.rld,
            rules.all.input.vsd,
            rules.all.input.csv,
            rules.all.input.heat,
            rules.all.input.pca,
            rules.all.input.vst,
            rules.all.input.rpl,
    log:    "LOGS/DE/run_deseq2.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD/2) if int(MAXTHREAD/2) >= 1 else 1
    params: bins   = BINS,
            outdir = lambda wildcards, output: os.path.dirname(outdir),
            compare = comparison
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DE/DESeq2_diffexp_2.R {input.anno} {input.cnt} {params.outdir} {params.compare} {threads} 2> {log} "

onsuccess:
    print("Workflow DE-deseq2 finished, no error")
