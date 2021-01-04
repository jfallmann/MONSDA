DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config2(SAMPLES,config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DE/DESEQ2/"
comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  plot = expand("{outdir}DESeq2_{comparison}_MA.pdf", outdir=outdir, comparison=compstr),
            tbl  = expand("{outdir}DESeq2_{comparison}.tsv.gz", outdir=outdir, comparison=compstr),
            sigtbl  = expand("{outdir}Sig_DESeq2_{comparison}.tsv.gz", outdir=outdir, comparison=compstr),
            heat = expand("{outdir}DESeq2_heatmap{i}.pdf", outdir=outdir,i=[1,2,3,"_samplebysample"]),
            pca  = expand("{outdir}DESeq2_PCA.pdf", outdir=outdir),
            vst  = expand("{outdir}DESeq2_VST_and_log2.pdf", outdir=outdir),
            rld  = expand("{outdir}DESeq2_rld.txt.gz", outdir=outdir),
            vsd  = expand("{outdir}DESeq2_vsd.txt.gz", outdir=outdir),
            session = expand("{outdir}DESeq2_SESSION.gz", outdir=outdir)# R object?

rule featurecount_unique:
    input:  reads = "MAPPED/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DE_deseq/{{file}}_tmp.counts", outdir=outdir)),
            cts   = "DE/Featurecounts_DE/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecounts_deseq2_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DE", DEENV)['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DE'),
             bins = BINS,
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_deseq2:
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: plt = rules.themall.input.plot,
            rld = rules.themall.input.rld,
            vsd = rules.themall.input.vsd,
            tbl = rules.themall.input.tbl,
            heat = rules.themall.input.heat,
            pca = rules.themall.input.pca,
            vst = rules.themall.input.vst,
            session = rules.themall.input.session
    log:    expand("LOGS/{outdir}run_deseq2.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEBIN]),
            outdir = outdir,
            compare = comparison
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {params.outdir} {params.compare} {threads} 2> {log}"

rule filter_significant_deseq2:
    input:  tbl = rules.run_deseq2.output.tbl
    output: sigtbl  = rules.themall.input.sigtbl
    log:    expand("LOGS/{outdir}filter_deseq2.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    shell: "set +o pipefail;for i in {outdir}DESeq2_*.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| head -n1 |gzip > {outdir}Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] <= -1.5 ||$F[2] >= 1.5) ){{print}}' |gzip >> {outdir}Sig_$fn && zcat $i| head -n1 |gzip > {outdir}SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] >= 1.5) ){{print}}' |gzip >> {outdir}SigUP_$fn && zcat $i|head -n1 |gzip > {outdir}SigDOWN_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] <= -1.5) ){{print}}' |gzip >> {outdir}SigDOWN_$fn; else touch {outdir}Sig_$fn {outdir}SigUP_$fn {outdir}SigDOWN_$fn; fi;done 2> {log}"
