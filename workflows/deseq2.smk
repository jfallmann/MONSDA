DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir = "DE/DESEQ2/"
comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  plot = expand("{outdir}{comparison}_DESeq2_MA.pdf", outdir=outdir, comparison=compstr),
            csv  = expand("{outdir}{comparison}_DESeq2.csv.gz", outdir=outdir, comparison=compstr),
            heat = expand("{outdir}DESeq2_heatmap{i}.pdf", outdir=outdir,i=[1,2,3,"_samplebysample"]),
            pca  = expand("{outdir}DESeq2_PCA.pdf", outdir=outdir),
            vst  = expand("{outdir}DESeq2_VST_and_log2.pdf", outdir=outdir),
            rld  = expand("{outdir}DESeq2_rld.txt.gz", outdir=outdir),
            vsd  = expand("{outdir}DESeq2_vsd.txt.gz", outdir=outdir),
            session = expand("{outdir}DESeq2_SESSION.gz", outdir=outdir)# R object?

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DE_deseq/{{file}}_tmp.counts", outdir=outdir)),
            cts   = expand("{outdir}Featurecounts_DE_deseq/{{file}}_mapped_sorted_unique.counts", outdir=outdir)
    log:    "LOGS/{file}/featurecounts_deseq2_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DE')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DE")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "snakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DE'),
             bins = BINS,
             tpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(samplecond(SAMPLES,config)[0], None ,config, "DE")['OPTIONS'][1].items())
    shell: "{params.bins}/Analysis/build_count_table_simple.py {params.dereps} --table {output.tbl} --anno {output.anno} {params.tpara} --loglevel DEBUG 2> {log}"

rule run_deseq2:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.themall.input.plot,
            rules.themall.input.rld,
            rules.themall.input.vsd,
            rules.themall.input.csv,
            rules.themall.input.heat,
            rules.themall.input.pca,
            rules.themall.input.vst,
            rules.themall.input.session
    log:    expand("LOGS/{outdir}run_deseq2.log",outdir=outdir)
    conda:  "snakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {params.outdir} {params.compare} {threads} 2> {log}"
