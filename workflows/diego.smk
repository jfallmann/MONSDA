DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/DIEGO/"
compare_string= comparable_as_string2(config,'DAS').split(",")
comparison=[i.split(":")[0] for i in compare_string]

rule themall:
    input:  dendrogram = expand("{outdir}{comparison}_dendrogram", outdir=outdir, comparison=comparison)

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

rule create_samplemaps:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output: smap = expand("{outdir}Tables/samplemap.txt", outdir=outdir),
            cmap = expand("{outdir}Tables/groupings.txt", outdir=outdir)
    log:    expand("LOGS/{outdir}create_samplemaps.log", outdir=outdir)
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd,config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd,config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"

rule prepare_junction_usage_matrix:
    input:  smap = rules.create_samplemaps.output.smap,
            cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output: tbl = expand("{outdir}Tables/junction_table_dexdas.txt.gz", outdir=outdir),
            anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:    expand("LOGS/{outdir}prepare_junction_usage_matrix.log", outdir=outdir)
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            dereps = lambda wildcards, input: get_reps(input.cnd,config,'DAS')
    shell:  "{params.bins}/Analysis/DAS/FeatureCounts2DIEGO.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule create_contrast_files:
    input:  cmap = rules.create_samplemaps.output.cmap
    output: contrast = expand("{outdir}Tables/{comparison}_contrast.txt", outdir=outdir, comparison=comparison)
    log:    expand("LOGS/{outdir}create_contrast_files.log", outdir=outdir)
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            compare=compare_string,
            outdir=outdir+'Tables/'
    shell:  "{params.bins}/Analysis/DAS/diego_contrast_files.py -g {input.cmap} -c {params.compare} -o {params.outdir} 2> {log}"

rule run_diego:
    input:  tbl = rules.prepare_junction_usage_matrix.output.tbl,
            contrast = expand(rules.create_contrast_files.output.contrast, outdir=outdir, comparison=comparison),
            # group = rules.create_samplemaps.output.cmap
    output: dendrogram = rules.themall.input.dendrogram,
            # grouplist = temp(expand("{outdir}subgroup", outdir=outdir))
    log:    expand("LOGS/{outdir}run_diego.log", outdir=outdir)
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "array1=({input.contrast}); array2=({output.dendrogram}); for i in ${{array[*]}}; do {params.bins} -a {input.tbl} -b ${{array1[$i]}} -x < (head -n 1 $i | awk '{print ${{array2[$i]}}}') -e -f {output.dendrogram} 2> {log}"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
