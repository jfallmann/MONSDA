DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/EDGER/"
comparison=comparable_as_string2(config,'DAS')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  all = expand("{outdir}All_Conditions_MDS.pdf", outdir=outdir),
            allsum = expand("{outdir}All_Conditions_sum_MDS.pdf", outdir=outdir),
            tbl = expand("{outdir}All_Conditions_normalized_table.tsv", outdir=outdir),
            bcv = expand("{outdir}All_Conditions_BCV.pdf", outdir=outdir),
            qld = expand("{outdir}All_Conditions_QLDisp.pdf", outdir=outdir),
            dift = expand("{outdir}{comparison}_diffSplice_{test}.tsv", outdir=outdir, comparison=compstr, test=["geneTest","simesTest","exonTest"]),
            tops = expand("{outdir}{comparison}_topSplice_simes_{n}.pdf", outdir=outdir, comparison=compstr, n=[str(i) for i in range(1,11)]),
            session = expand("{outdir}EDGER_DAS_SESSION.gz", outdir=outdir)

rule featurecount_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DEU_edger/{{file}}_tmp.counts", outdir=outdir)),
            cts   = expand("{outdir}Featurecounts_DEU_edger/{{file}}_mapped_sorted_unique.counts", outdir=outdir)
    log:    "LOGS/{file}/featurecount_DAS_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count = COUNTBIN,
            anno  = lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']]),
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DAS")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DAS'),
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
            rules.themall.input.tops,
            rules.themall.input.session
    log:    expand("LOGS/{outdir}run_edger.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DASBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.outdir} {params.compare} {threads} 2> {log} "
