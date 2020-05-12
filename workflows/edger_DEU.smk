DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DEU/EDGER/"
comparison=comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  all = expand("{outdir}All_Conditions_MDS.pdf", outdir=outdir),
            allsum = expand("{outdir}All_Conditions_sum_MDS.pdf", outdir=outdir),
            tbl = expand("{outdir}All_Conditions_normalized_table.tsv", outdir=outdir),
            bcv = expand("{outdir}All_Conditions_BCV.pdf", outdir=outdir),
            qld = expand("{outdir}All_Conditions_QLDisp.pdf", outdir=outdir),
            dift = expand("{outdir}{comparison}_exons_{sort}.tsv", outdir=outdir, comparison=compstr, sort=["logFC-sorted","pValue-sorted"]),
            plot = expand("{outdir}{comparison}_MD.pdf", outdir=outdir, comparison=compstr),
            session = expand("{outdir}EDGER_DEU_SESSION.gz", outdir=outdir)

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: cts   = expand("{outdir}Featurecounts_DEU_edger/{{file}}_mapped_sorted_unique.counts", outdir=outdir)
    log:    "LOGS/{file}/featurecounts_DEU_edger_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DEU')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.cts} {input.reads} 2> {log}"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS
    shell: "{params.bins}/Analysis/build_count_table_id.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edger:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: rules.themall.input.all,
            rules.themall.input.allsum,
            rules.themall.input.tbl,
            rules.themall.input.bcv,
            rules.themall.input.qld,
            rules.themall.input.dift,
            rules.themall.input.plot,
            rules.themall.input.session
    log:    expand("LOGS/{outdir}run_edger.log",outdir=outdir)
    conda:  "snakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.outdir} {params.compare} {threads} 2> {log}"
