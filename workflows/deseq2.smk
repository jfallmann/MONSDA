DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,'DE')
COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

rule all:
    input:  "DE/DESEQ2/DONE"

rule featurecount_unique:
    input:  reads = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: tbl = "COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts",
            anno = temp("COUNTS/Featurecounter_genes/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_de_gene_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t gene -g '+config['COUNTING']['FEATURES']['gene'],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output.anno} -o {output.tbl} {input.reads} 2> {log}"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_unique.output.tbl, file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/RUN_DE_Analysis.tbl.gz",
             anno = "DE/Tables/RUN_DE_Analysis.anno.gz"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  decond = lambda wildcards, input: str.join(',',[','.join(tool_params_rep(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITIONS']) for x in [str.join(os.sep,y.split(os.sep)[:-1]) for y in input.cnd]]),
             dereps = lambda wildcards, input: str.join(',',[str.join(os.sep, str.join(os.sep,x.split(os.sep)[2:-1])+',',z) for z in [tool_params_rep(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['REPLICATES'] for x in [str.join(os.sep,y.split(os.sep)[:-1]) for y in input.cnd]]]),
             detypes = lambda wildcards, input: '-t '+str.join(',',[','.join(tool_params_rep(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['TYPES']) for x in [str.join(os.sep,y.split(os.sep)[:-1]) for y in input.cnd]]),
             #samples = lambda wildcards, input: str.join(',',list(set(input.cnd))),
             paired = lambda wildcards, input:  str.join(',',[checkpaired_rep([str.join(os.sep,x.split(os.sep)[2:])],config) for x in list(set(input.cnd))]),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py -r {params.dereps} -c {params.decond} {params.detypes} --paired {params.paired} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_deseq2:
    input:  cnt = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: csv = "DE/DESEQ2/DONE"
    log:    "LOGS/DE/run_deseq2.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD/2) if int(MAXTHREAD/2) >= 1 else 1
    params: bins = BINS,
            outdir = lambda wildcards, output: os.path.dirname(output.csv),
            #condcombs = lambda wildcards, input: ','.join([map(str, comb) for comb in combinations([','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITION']) for x in input.cnt],2)]),
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DE/DESeq2_diffexp.R {input.anno} {input.cnt} {params.outdir} {threads} 2> {log} && touch {output.csv}"

#rule themall:
#    input:  rules.summarize_counts.output
#    output: "COUNTS/DONE"
#    conda:  "snakes/envs/base.yaml"
#    threads: 1
#    params: bins = BINS
#    shell:  "touch {output}"

onsuccess:
    print("Workflow finished, no error")
