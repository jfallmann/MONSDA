DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,'DE')

rule all:
    input:  "DE/DESEQ2/DONE"

rule prepare_count_table:
    input:   cnd = expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/RUN_DE_Analysis.tbl",
             anno = "DE/Tables/RUN_DE_Analysis.anno"
    log:     "LOGS/DE/prepare_count_table.log"
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
    output: csv = "DE/DESEQ2/DONE"
    log:    "LOGS/DE/run_deseq2.log"
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
