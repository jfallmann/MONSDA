for analysis in ['DE', 'DEU', 'DAS']:
    if analysis in config:
        DEUBIN, DEUENV = env_bin_from_config2(SAMPLES,config,analysis)
        COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')
    else:
        continue

    rule all:
        input:  analysis+"/DEXSEQ/DONE"

    if analysis == 'DEU':
        rule featurecount_unique:
            input:  mapf = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
            output: cts = "COUNTS/Featurecounter_dexseq/{file}_mapped_sorted_unique.counts",
                    anno = "COUNTS/Featurecounter_dexseq/{file}_dexseq.gff.gz"
                    #tan = temp("COUNTS/Featurecounter_dexseq/{file}_anno.tmp")
            log:    "LOGS/{file}/featurecount_"+analysis+"_dexseq_unique.log"
            conda:  "snakes/envs/"+COUNTENV+".yaml"
            threads: MAXTHREAD
            params: count = COUNTBIN,
                    anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]),
                    cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items()),
                    paired = lambda x: '-p' if paired == 'paired' else '',
                    stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
                    bins = BINS,
                    countgtf = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_fc_dexseq.gtf')),
                    dexgtf = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_dexseq.gtf'))
            shell:  "if [ ! -f \"{params.dexgtf}\"] || [ ! -f \"{params.countgtf}\" ];then {params.bins}/Analysis/DEU/prepare_dexseq_annotation2.py -f {params.countgtf} {params.anno} {params.dexgtf} ;fi && ln -s {params.dexgtf} {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.countgtf}) -o {output.cts} {input.mapf} 2> {log}"

        rule prepare_count_table:
            input:   cnd = expand("COUNTS/Featurecounter_dexseq/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
            output:  tbl = "DEU/Tables/RUN_DEU_Analysis.tbl.gz",
                     anno = "DEU/Tables/RUN_DEU_Analysis.anno.gz"
            log:     "LOGS/DEU/prepare_count_table.log"
            conda:   "snakes/envs/"+DEUENV+".yaml"
            threads: 1
            params:  decond = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DEU')['CONDITION']) for x in input.cnd]),
                     dereps = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DEU')['REPLICATES']) for x in input.cnd]),
                     samples = lambda wildcards, input: str.join(',',input.cnd),
                     bins = BINS
            shell: "{params.bins}/Analysis/DEU/build_DEXSeq_table.py -l {params.samples} -r {params.dereps} -c {params.decond} --table {output.tbl} --anno {output.anno} 2> {log}"

        rule run_dexeq:
            input:  cnt = rules.prepare_count_table.output.tbl,
                    anno = rules.prepare_count_table.output.anno,
                    flat = expand(rules.featurecount_unique.output.anno,file=samplecond(SAMPLES,config))
            output: csv = "DEU/DEXSEQ/DONE"
            log:    "LOGS/DEU/run_deseq2.log"
            conda:  "snakes/envs/"+DEUENV+".yaml"
            threads: 1
            params: bins = BINS,
                    outdir = lambda wildcards, output: os.path.dirname(output.csv),
                    #condcombs = lambda wildcards, input: ','.join([map(str, comb) for comb in combinations([','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['CONDITION']) for x in input.cnt],2)]),
            shell: "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DEU/DEXSeq.R {input.anno} {input.cnt} {params.outdir} 2> {log} && touch {output.csv}"

#rule themall:
#    input:  rules.summarize_counts.output
#    output: "COUNTS/DONE"
#    conda:  "snakes/envs/base.yaml"
#    threads: 1
#    params: bins = BINS
#    shell:  "touch {output}"

    onsuccess:
        print("Workflow "+analysis+"finished, no error")
