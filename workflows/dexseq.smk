DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DEU/DEXSEQ/"
comparison = comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input: tbl  = expand("{outdir}DEXSeq_{comparison}.tsv.gz", outdir=outdir, comparison=compstr),
           sigtbl  = expand("{outdir}Sig_DEXSeq_{comparison}.tsv.gz", outdir=outdir, comparison=compstr),
           plot = expand("{outdir}DEXSeq_{comparison}_DispEsts.pdf", outdir=outdir, comparison=compstr),
           html = expand("{outdir}DEXSeqReport_{comparison}/DEXSeq_{comparison}.html", outdir=outdir, comparison=compstr),
           session = expand("{outdir}DEXSeq_SESSION.gz", outdir=outdir)

rule prepare_deu_annotation:
    input:   anno = ANNOTATION
    output:  countgtf = expand("{refd}/{countanno}", ref=REFDIR, countanno=ANNOTATION.replace('.gtf','_fc_dexseq.gtf')),
             deugtf   = expand("{refd}/{deuanno}", ref=REFDIR, deuanno=ANNOTATION.replace('.gtf','_dexseq.gtf'))
    log:     "LOGS/featurecount_dexseq_annotation.log"
    conda:   "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params:  bins = BINS,
             countstrand = lambda x: '-s' if stranded == 'fr' or stranded == 'rf' else ''
    shell:  "{params.bins}/Analysis/DEU/prepare_deu_annotation2.py -f {output.countgtf} {params.countstrand} {input.anno} {output.deugtf} 2>> {log}"

rule featurecount_dexseq_unique:
    input:  reads = "MAPPED/{combo}{file}_mapped_sorted_unique.bam",
            countgtf = expand(rules.prepare_deu_annotation.output.countgtf, ref=REFDIR, countanno=ANNOTATION.replace('.gtf','_fc_dexseq.gtf')),
            deugtf = expand(rules.prepare_deu_annotation.output.deugtf, ref=REFDIR, deuanno=ANNOTATION.replace('.gtf','_dexseq.gtf'))
    output: tmp   = temp(expand("{outdir}Featurecounts_DEU_dexseq/{{file}}_tmp.counts", outdir=outdir)),
            cts   = "DEU/Featurecounts_DEU_dexseq/{combo}{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{combo}{file}/featurecounts_dexseq_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb  = COUNTBIN,
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU", COUNTENV)['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {input.countgtf}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k2,2 -k3,3n -k4,4n -k1,1 -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_dexseq_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS,
    shell: "{params.bins}/Analysis/DEU/build_DEU_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_dexseq:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
            flat = rules.prepare_deu_annotation.output.deugtf
    output: plot = rules.themall.input.plot,
            tbl  = rules.themall.input.tbl,
            html = rules.themall.input.html,
            session = rules.themall.input.session
    log:    expand("LOGS/{outdir}run_dexseq.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {input.flat} {params.outdir} {params.compare} {threads} 2> {log}"

rule filter_significant_dexseq:
    input:  tbl = rules.run_dexseq.output.tbl
    output: sigtbl  = rules.themall.input.sigtbl
    log:    expand("LOGS/{outdir}filter_dexseq.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    shell: "set +o pipefail;for i in {outdir}DEXSeq_*.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i|head -n1 |gzip > {outdir}Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < 0.05 && ($F[9] <= -1.5 ||$F[9] >= 1.5) ){{print}}' |gzip >> {outdir}Sig_$fn && zcat $i|head -n1 |gzip > {outdir}SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < 0.05 && ($F[9] >= 1.5) ){{print}}' |gzip >> {outdir}SigUP_$fn && zcat $i|head -n1 |gzip > {outdir}SigDOWN_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < 0.05 && ($F[9] <= -1.5) ){{print}}' |gzip >> {outdir}SigDOWN_$fn;else touch {outdir}Sig_$fn {outdir}SigUP_$fn {outdir}SigDOWN_$fn; fi;done 2> {log}"
