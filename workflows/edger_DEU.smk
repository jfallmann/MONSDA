DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config2(SAMPLES,config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DEU/EDGER/"
comparison = comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  all = expand("{outdir}EDGER_DEU_All_Conditions_MDS.png", outdir=outdir),
            tbl = expand("{outdir}EDGER_DEU_All_Conditions_normalized.tsv.gz", outdir=outdir),
            bcv = expand("{outdir}EDGER_DEU_All_Conditions_BCV.png", outdir=outdir),
            qld = expand("{outdir}EDGER_DEU_All_Conditions_QLDisp.png", outdir=outdir),
            dift = expand("{outdir}EDGER_DEU_{comparison}_exons_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["logFC-sorted","pValue-sorted"]),
            sigdift = expand("{outdir}Sig_EDGER_DEU_{comparison}_exons_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["pValue-sorted"]),
            plot = expand("{outdir}EDGER_DEU_{comparison}_MD.png", outdir=outdir, comparison=compstr),
            session = expand("{outdir}EDGER_DEU_SESSION.gz", outdir=outdir)

rule featurecount_unique:
    input:  reads = "MAPPED/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DEU_edger/{{file}}_tmp.counts", outdir=outdir)),
            cts   = "DEU/Featurecounts_DEU/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{file}/featurecounts_DEU_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno  = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEU")['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output:  tbl  = expand("{outdir}Tables/COUNTS.gz",outdir=outdir),
             anno = expand("{outdir}Tables/ANNOTATION.gz",outdir=outdir)
    log:     expand("LOGS/{outdir}prepare_count_table.log",outdir=outdir)
    conda:   "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd,config,'DEU'),
             bins = BINS,
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edgerDEU:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: all = rules.themall.input.all,
            tbl = rules.themall.input.tbl,
            bcv = rules.themall.input.bcv,
            qld = rules.themall.input.qld,
            dift = rules.themall.input.dift,
            plot = rules.themall.input.plot,
            session = rules.themall.input.session
    log:    expand("LOGS/{outdir}run_edger.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEUBIN]),
            outdir = outdir,
            compare = comparison
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.outdir} {params.compare} {threads} 2> {log}"

rule filter_significant_edgerDEU:
    input:  dift = rules.run_edgerDEU.output.dift
    output: sigdift = rules.themall.input.sigdift
    log:    expand("LOGS/{outdir}filter_edgerDEU.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    shell: "set +o pipefail;for i in {outdir}EDGER_DEU*pValue-sorted.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i|head -n1 |gzip > {outdir}Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] <= -1.5 ||$F[2] >= 1.5) ){{print}}' |gzip >> {outdir}Sig_$fn && zcat $i|head -n1 |gzip > {outdir}SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] >= 1.5) ){{print}}' |gzip >> {outdir}SigUP_$fn && zcat $i|head -n1 |gzip > {outdir}SigDOWN_$fn &&zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[6]);if ($F[6] < 0.05 && ($F[2] <= -1.5) ){{print}}' |gzip >> {outdir}SigDOWN_$fn; else touch {outdir}Sig_$fn {outdir}SigUP_$fn {outdir}SigDOWN_$fn; fi;done 2> {log}"
