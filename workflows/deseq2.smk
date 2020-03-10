DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,'DE')
COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

rule all:
    input:  "DE/DESEQ2/DONE"

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = "COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecount_de_gene_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t gene -g '+config['COUNTING']['FEATURES']['gene'],
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = "DE/Tables/RUN_DE_Analysis.tbl.gz",
             anno = "DE/Tables/RUN_DE_Analysis.anno.gz"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DE'),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_deseq2:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: csv  = "DE/DESEQ2/DONE"
    log:    "LOGS/DE/run_deseq2.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD/2) if int(MAXTHREAD/2) >= 1 else 1
    params: bins   = BINS,
            outdir = lambda wildcards, output: os.path.dirname(output.csv),
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DE/DESeq2_diffexp.R {input.anno} {input.cnt} {params.outdir} {threads} 2> {log} && touch {output.csv}"

onsuccess:
    print("Workflow DE finished, no error")
