for analysis in ['DE', 'DEU', 'DAS']:
    if analysis in config:
        DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,analysis)
        COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')


    rule all:
        input:  analysis+"/EDGER/DONE"

    if analysis == 'DE':

        rule featurecount_unique:
            input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: "COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts",
                    temp("COUNTS/Featurecounter_genes/{file}_unique.anno")
            log:    "LOGS/{file}/featurecount_de_gene_unique.log"
            conda:  "snakes/envs/"+COUNTENV+".yaml"
            threads: MAXTHREAD
            params: count = COUNTBIN,
                    anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
                    cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items()),
                    paired = lambda x: '-p' if paired == 'paired' else '',
                    stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
                    outfile = lambda wildcards: expand("COUNTS/Featurecounter_{region}/{file}_mapped_sorted.counts",file = wildcards.file, region = 'gene')
            shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output[1]} -o {params.outfile} {input[0]} 2> {log}"

        rule prepare_count_table:
            input:   cnd = expand("COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
            output:  tbl = "DE/EDGER/Tables/RUN_DE_Analysis.tbl.gz",
                     anno = "DE/EDGER/Tables/RUN_DE_Analysis.anno.gz"
            log:     "LOGS/DE/EDGER/prepare_count_table.log"
            conda:   "snakes/envs/"+DEENV+".yaml"
            threads: 1
            params:  decond = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITION']) for x in input.cnd]),
                     dereps = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['REPLICATES']) for x in input.cnd]),
                     samples = lambda wildcards, input: str.join(',',input.cnd),
                     bins = BINS
            shell: "{params.bins}/Analysis/DE/build_DESeq_table.py -l {params.samples} -r {params.dereps} -c {params.decond} --table {output.tbl} --anno {output.anno} 2> {log}"

        rule run_deseq2:
            input:  cnt = rules.prepare_count_table.output.tbl,
                    anno = rules.prepare_count_table.output.anno,
            output: csv = "DE/EDGER/DONE"
            log:    "LOGS/DE/run_edger.log"
            conda:  "snakes/envs/"+DEENV+".yaml"
            threads: 1
            params: bins = BINS,
                    outdir = lambda wildcards, output: os.path.dirname(output.csv),
                    #condcombs = lambda wildcards, input: ','.join([map(str, comb) for comb in combinations([','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITION']) for x in input.cnt],2)]),
            shell: "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DE/DESeq2_diffexp.R {input.anno} {input.cnt} {params.outdir} 2> {log} && touch {output.csv}"

#rule themall:
#    input:  rules.summarize_counts.output
#    output: "COUNTS/DONE"
#    conda:  "snakes/envs/base.yaml"
#    threads: 1
#    params: bins = BINS
#    shell:  "touch {output}"

onsuccess:
    print("Workflow finished, no error")
