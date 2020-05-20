DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DEU/DEXSEQ/"
comparison=comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input: tbl  = expand("{outdir}DEXSeq_{comparison}.tsv.gz", outdir=outdir, comparison=compstr),
           plot = expand("{outdir}DEXSeq_{comparison}_DispEsts.pdf", outdir=outdir, comparison=compstr),
           html = expand("{outdir}DEXSeqReport_{comparison}/DEXSeq_{comparison}.html", outdir=outdir, comparison=compstr),
           session = expand("{outdir}DEXSeq_SESSION.gz", outdir=outdir)

rule prepare_count_annotation:
    input:   anno   = expand("{ref}/{gen}/{anno}", ref=REFERENCE, gen=os.path.dirname(genomepath(SAMPLES[0],config)), anno=tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION'])
    output:  countgtf = expand("{ref}/{gen}/{countanno}", ref=REFERENCE, gen=os.path.dirname(genomepath(SAMPLES[0],config)), countanno=tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION'].replace('.gtf','_fc_dexseq.gtf')),
             deugtf   = expand("{ref}/{gen}/{deuanno}", ref=REFERENCE, gen=os.path.dirname(genomepath(SAMPLES[0],config)), deuanno=tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION'].replace('.gtf','_dexseq.gtf'))
    log:     "LOGS/featurecount_dexseq_annotation.log"
    conda:   "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params:  bins = BINS,
             countstrand = lambda x: '-s' if stranded == 'fr' or stranded == 'rf' else ''
    shell:  "{params.bins}/Analysis/DEU/prepare_deu_annotation2.py -f {output.countgtf} {params.countstrand} {input.anno} {output.deugtf} 2>> {log}"

rule featurecount_dexseq_unique:
    input:  reads = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            countgtf = expand(rules.prepare_count_annotation.output.countgtf, ref=REFERENCE, gen=os.path.dirname(genomepath(SAMPLES[0],config)), countanno=tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION'].replace('.gtf','_fc_dexseq.gtf')),
            deugtf = expand(rules.prepare_count_annotation.output.deugtf, ref=REFERENCE, gen=os.path.dirname(genomepath(SAMPLES[0],config)), deuanno=tool_params(SAMPLES[0], None, config, 'DEU')['ANNOTATION'].replace('.gtf','_dexseq.gtf'))
    output: tmp   = temp(expand("{outdir}Featurecounts_DEU_edger/{{file}}_tmp.counts", outdir=outdir)),
            cts   = expand("{outdir}Featurecounts_DEU_edger/{{file}}_mapped_sorted_unique.counts", outdir=outdir)
    log:    "LOGS/{file}/featurecounts_dexseq_unique.log"
    conda:  "snakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: count  = COUNTBIN,
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU")['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.count} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {input.countgtf}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_dexseq_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = expand("{outdir}Tables/RUN_DEU_Analysis.tbl.gz",outdir=outdir),
             anno = expand("{outdir}Tables/RUN_DEU_Analysis.anno.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "snakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS
    shell: "{params.bins}/Analysis/DEU/build_DEU_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_dexseq:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
            flat = rules.prepare_count_annotation.output.deugtf
    output: plot = rules.themall.input.plot,
            tbl  = rules.themall.input.tbl,
            html = rules.themall.input.html,
            session = rules.themall.input.session
    log:    expand("LOGS/{outdir}run_dexseq.log",outdir=outdir)
    conda:  "snakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {input.flat} {params.outdir} {params.compare} {threads} 2> {log}"
