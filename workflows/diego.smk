DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/DIEGO/"
comparison=comparable_as_string2(config,'DAS')
comps = comparison.split(",")
rule themall:
    input:  expand("{outdir}dendrogram", outdir=outdir)

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = "COUNTS/Featurecounts_DAS_diego/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecounts_DAS_diego_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DAS")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule create_samplemap:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output: smap = expand("{outdir}Tables/samplemap.txt",outdir=outdir),
            cmap = expand("{outdir}Tables/groupings.txt",outdir=outdir)
    log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd,config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd,config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"

rule prepare_junction_usage_matrix:
    input:  smap  = rules.create_samplemaps.output.smap
    output: tbl  = expand("{outdir}Tables/junction_table_dexdas.txt",outdir=outdir)
    log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS
    shell:  "perl {params.bins}/Analysis/DAS/FeatureCounts2DIEGO.pl -i {input.smap} -o {output.tbl} 2> {log}"

rule run_diego:
    input:  tbl= rules.prepare_junction_usage_matrix.output.tbl,
            group = rules.create_samplemap.output.cmap
    output: expand("{outdir}dendrogram", outdir=outdir)
    log:    "LOGS/"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "python {params.bins} -a {input.tbl} -b {input.group} -x <(head -n 1 {input.group} | awk '{print $1}') -e -f {output}"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
