
DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DE/EDGER/"
comparison=comparable_as_string(config,'DE')

rule all:
    input:  plot = expand("{outdir}{comparison}.png", outdir=outdir, comparison=comparison.split(",")),
            bcv = expand("{outdir}BCV.png", outdir=outdir),
            mds = expand("{outdir}MDS.png", outdir=outdir),
            tbl = expand("{outdir}normalized_table.tsv", outdir=outdir)

rule featurecount_unique:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: "COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts",
            temp("COUNTS/Featurecounter_genes/{file}_unique.anno")
    log:    "LOGS/{file}/featurecount_de_gene_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'COUNTING')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "COUNTING")['OPTIONS'][0].items())+' -t gene -g '+config['COUNTING']['FEATURES']['gene'],
            paired = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "zcat {params.anno} > {output[1]} && {params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a {output[1]} -o {output[0]} {input[0]} 2> {log}"

rule prepare_count_table:
    input:   cnd = expand("COUNTS/Featurecounter_genes/{file}_mapped_sorted_unique.counts", file=samplecond(SAMPLES,config))
    output:  tbl = "DE/Tables/RUN_DE_Analysis.tbl.gz",
             anno = "DE/Tables/RUN_DE_Analysis.anno.gz"
    log:     "LOGS/DE/prepare_count_table.log"
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  decond = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')["GROUP"]) for x in input.cnd]),
             dereps = lambda wildcards, input: str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['REPLICATES']) for x in input.cnd]),
             detypes = lambda wildcards, input: '-t '+str.join(',',[','.join(tool_params(str.join(os.sep, x.split(os.sep)[2:]).replace('_mapped_sorted_unique.counts',''), None, config, 'DE')['TYPES']) for x in input.cnd]),
             samples = lambda wildcards, input: str.join(',',input.cnd),
             bins = BINS
    shell: "{params.bins}/Analysis/DE/build_DESeq_table.py -l {params.samples} -r {params.dereps} -c {params.decond} -t {params.detypes} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_edger:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.all.input.plot,
            rules.all.input.bcv,
            rules.all.input.mds,
            rules.all.input.tbl
    log:    "LOGS/DE/run_edger.log"
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params: bins = BINS,
            outdir = lambda wildcards, output: os.path.dirname(outdir),
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins}/Analysis/DE/EdgeR.R {input.tbl} {input.anno} {params.outdir} {params.compare} 2> {log} "


onsuccess:
    print("Workflow DE-edgeR finished, no error")
