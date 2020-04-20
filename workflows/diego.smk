DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/DIEGO/"
comparison=comparable_as_string2(config,'DAS')

rule themall:
    input:  expand("{outdir}{comparison}_dendrogram", outdir=outdir, comparison=[i.split(":")[0] for i in comparison.split(",")])

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = expand("COUNTS/Featurecounts_DAS_diego/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
    log:    "LOGS/{file}/featurecounts_DAS_diego_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DAS")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule create_samplemaps:
    input:  cnd  = rules.featurecount_unique.output.cts
    output: smap = "{outdir}Tables/samplemap.txt",
            cmap = "{outdir}Tables/groupings.txt"
    log:    "LOGS/DAS/DIEGO/create_samplemaps.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd,config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd,config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"

rule prepare_junction_usage_matrix:
    input:  smap = rules.create_samplemaps.output.smap
    output: tbl = "{outdir}Tables/junction_table_dexdas.txt"
    log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            dereps = lambda wildcards, input: get_reps(input.cnd,config,'DAS')
    shell:  "{params.bins}/Analysis/DAS/FeatureCounts2DIEGO.py {params.dereps} --table {output.tbl}  --anno {output.anno} 2> {log}"

rule create_contrast_files:
    input:  rules.create_samplemaps.output.cmap
    output: "{outdir}{comparison}_contrast.txt"
    log:    "LOGS/DAS/DIEGO/create_contrast_files.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            compare=comparison,
            outdir=outdir
    shell:  "python3 {params.bins}/Analysis/DAS/diego_contrast_files.py -g {input} -c {params.compare} -o {params.outdir} 2> {log}"

rule run_diego:
    input:  tbl = rules.prepare_junction_usage_matrix.output.tbl,
            contrast = rules.create_contrast_files.output
    output: rules.themall.input,
            grouplist = temp(expand("{outdir}subroup", outdir=outdir))
    log:    "LOGS/DAS/DIEGO/run_diego.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "head -n 1 {input.group} | awk '{{print $1}}' > {output.grouplist} && python {params.bins} -a {input.tbl} -b {input.contrast} -x {output.grouplist} -e -f {output} 2> {log}"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
