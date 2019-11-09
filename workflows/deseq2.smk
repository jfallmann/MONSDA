DEBIN, DEENV = env_bin_from_config2(SAMPLES,config,'DE')
combilist = [os.path.join(x) for x in tool_params(input[0], None, config, 'DE')['COMBINATIONS'])]  # Need to adjust for list of lists

#get files with specified pattern
one = [os.path.abspath(os.path.join(x, "*_mapped_sorted_unique.counts")) for x in combilist[0]]
two = [os.path.abspath(os.path.join(x, "*_mapped_sorted_unique.counts")) for x in combilist[1]]

#search for files
cond_one = [re.sub('\_mapped\_sorted\_unique\.counts', '', x) for x in natsorted(glob.glob(one), key=lambda y: y.lower())]
cond_two = [re.sub('\_mapped\_sorted\_unique\.counts', '', x) for x in natsorted(glob.glob(two), key=lambda y: y.lower())]

rule all:
    input:  expand("DE/DESEQ2/{file}.csv",file=samplecond(SAMPLES,config)),

rule prepare_count_table:
    input:   cndo = expand({file}_mapped_sorted_unique.counts, file=cond_one),
             cndt = expand({file}_mapped_sorted_unique.counts, file=cond_two)
    output:  tbl = "DE/Tables/{file}.tbl",
             tc  = temp("DE/Tables/{file}.csv")
    log:     "LOGS/DE/{file}.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  decombi = lambda wildcards, input: tool_params(input[0], None, config, 'DE')['COMBINATIONS']),
             bins = BINS,
             co = combilist[0]
             ct = combilist[1]
    shell: "for i in {input.cndo}; do echo \"{params.co} cond_one \$i\"> {output.tc};done && for i in {input.cndt}; do echo \"{params.ct} cond_one \$i\">> {output.tc};done && python2 {params.bins}/Analysis/DE/build_DESeq_table.py -l {output.tc} -n > {output.tbl} 2> {log}"

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
