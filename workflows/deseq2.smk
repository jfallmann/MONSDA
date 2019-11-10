DEBIN, DEENV = list(env_bin_from_config2(SAMPLES,config,'DE'))[0:2]

rule all:
    input:  expand("DE/DESEQ2/{file}.csv",file=samplecond(SAMPLES,config)),

rule prepare_count_table:
    input:   cnd = expand("COUNTS/Featurecounter/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/RUN_DE_Analysis.tbl"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  decond = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITION']) for x in input.cnd]),
             dereps = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['REPLICATES']) for x in input.cnd]),
             samples = lambda wildcards, input: str.join(',',input.cnd),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py -l {params.samples} -r {params.dereps} -c {params.decond} 1> {output.tbl} 2> {log}"

rule run_deseq2:
    input:  cnt = expand(rules.prepare_count_table.output, file=samplecond(SAMPLES,config))
    output: csv = "DE/DESEQ2/{file}.csv"
    log:    "LOGS/DE/{file}.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params: bins = BINS
    shell: "Rscript {params.bins}/Analysis/DE/DESeq2_diffexp.R {input.cnt} {output.csv} 2> {log}"

#rule themall:
#    input:  rules.summarize_counts.output
#    output: "COUNTS/DONE"
#    conda:  "snakes/envs/base.yaml"
#    threads: 1
#    params: bins = BINS
#    shell:  "touch {output}"

onsuccess:
    print("Workflow finished, no error")
