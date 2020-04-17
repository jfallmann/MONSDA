DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/DIEGO/"
comparison=comparable_as_string2(config,'DAS')
comps = comparison.split(",")


rule themall:
    input: tbl = expand("{outdir}DIEGO_{comparison}.tsv.gz", outdir=outdir, comparison=comparison.split(",")),
           plot = expand("{outdir}DIEGO_{comparison}_DispEsts.pdf", outdir=outdir, comparison=comparison.split(",")),
           html = expand("{outdir}DIEGO_{comparison}/DEXSeq_{comparison}.html", outdir=outdir, comparison=comparison.split(",")),
           session = expand("{outdir}DEXSeq_SESSION.gz", outdir=outdir)# R object?

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = "COUNTS/Featurecounter_DAS_diego/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecount_DAS_diego_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DAS")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule create_genome_annotation_file:
    input:  gff = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']])
    output: bed = expand("{outdir}Tables/Annotation_DIEGO.bed", outdir=outdir)
    log:    "LOGS/DAS/DIEGO/create_genome_annotation_file.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    shell:  "perl gfftoDIEGObed.pl -g  <(perl -F\\\\040 -wlane '{($F[0] = $F[0] =~ /^chr/ ? $F[0] : \"chr\".$F[0])=~ s/\_/\./g;print $F[0]}' <(zcat {input.gff})) -o {output.bed} 2> {log}"

rule create_samplemap:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output: smap  = expand("{outdir}Tables/samplemap.txt",outdir=outdir),
            cmap  = expand("{outdir}Tables/groupings.txt",outdir=outdir)
    log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd,config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd,config,'DAS'),
            bins = BINS
    shell:  "echo {params.slist} > {output.smap} && echo {params.clist} > {output.cmap} 2> {log}"

rule prepare_junction_usage_matrix:
    input:  smap  = rules.create_samplemap.output.smap
    output: tbl  = expand("{outdir}Tables/junction_table_dexdas.txt",outdir=outdir)
    log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: 1
    shell: " perl HTseq2DIEGO.pl -i {input.smap} -o {output.tbl} 2> {log}"

rule run_diego:
    input:  tbl= rules.prepare_junction_usage_matrix.output.tbl,
            anno = rules.create_genome_annotation_file.output.bed,
            base = rules.create_samplemap.output.cmap
    output: expand("{outdir}dendrogram", outdir=outdir)
    log:    "LOGS/"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "python {params.bins} -a {input.tbl} -b {input.anno} -x <(head -n 1 {input.base} | awk '{print $1}') -e -f {output}"

onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
