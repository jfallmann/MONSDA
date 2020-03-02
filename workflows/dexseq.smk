DEUBIN, DEUENV = env_bin_from_config2(SAMPLES,config,'DEU')
COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

rule all:
    input:  "DEU/DEXSEQ/DONE"

rule featurecount_unique:
    input:  mapf = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: cts  = "COUNTS/Featurecounter_dexseq/{file}_mapped_sorted_unique.counts",
            anno = "COUNTS/Featurecounter_dexseq/{file}_dexseq.gtf.gz"
    log:    "LOGS/{file}/featurecount_dexseq_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count  = COUNTBIN,
            anno   = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]),
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            countstrand = lambda x: 'True' if stranded == 'fr' or stranded == 'rf' else 'False',
            countgtf = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_fc_dexseq.gtf')),
            dexgtf   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_dexseq.gtf'))
    shell:  "if [ ! -f \"{params.dexgtf}\" ] || [ ! -f \"{params.countgtf}\" ];then {params.bins}/Analysis/DEU/prepare_dexseq_annotation2.py -f {params.countgtf} -s {params.countstrand} {params.anno} {params.dexgtf} ;fi && ln -s {params.dexgtf} {output.anno} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.countgtf}) -o {output.cts} {input.mapf} 2> {log}"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = "DEU/Tables/RUN_DEU_Analysis.tbl.gz",
             anno = "DEU/Tables/RUN_DEU_Analysis.anno.gz"
    log:     "LOGS/DEU/prepare_count_table.log"
    conda:   "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  decond  = lambda wildcards, input: str.join(',',get_reps(input.cnd,config,'DEU','CONDITIONS')),
             dereps  = lambda wildcards, input: str.join(',',get_reps(input.cnd,config,'DEU','REPLICATES')),
             detypes = lambda wildcards, input: '-t '+str.join(',',get_reps(input.cnd,config,'DEU','TYPES')),
             paired  = lambda wildcards, input:  '--paired '+str.join(',',[checkpaired_rep([str.join(os.sep,x.split(os.sep)[2:]) for x in get_reps(input.cnd,config,'DEU','REPLICATES')],config)]),
             bins = BINS
    shell: "{params.bins}/Analysis/DEU/build_DEXSeq_table.py -r {params.dereps} -c {params.decond} {params.detypes} {params.paired} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_dexeq:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno
    output: csv  = "DEU/DEXSEQ/DONE"
    log:    "LOGS/DEU/run_deseq2.log"
    conda:  "snakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD/2) if int(MAXTHREAD/2) >= 1 else 1
    params: bins   = BINS,
            outdir = lambda wildcards, output: os.path.dirname(output.csv),
            flat   = lambda wildcards: os.path.abspath(str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(SAMPLES[0], config)),tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION']]).replace('.gtf','_dexseq.gtf'))
    shell: "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DEU/DEXSeq.R {input.anno} {input.cnt} {params.flat} {params.outdir} {threads} 2> {log} && touch {output.csv}"

onsuccess:
    print("Workflow DEU finished, no error")
