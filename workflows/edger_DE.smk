DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

combi = "COMBINATION"

outdir = "DE/EDGER/"
comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("{outdir}/DE_EDGER_{combi}_SESSION.gz", outdir=outdir, combi=combi),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/SUM_DE_EDGER.Rmd")

rule featurecount_unique:
    input:  reads = "MAPPED/{combo}{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DE_edger/{{file}}_tmp.counts", outdir=outdir)),
            cts   = "DE/Featurecounts_DE/{combo}{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{combo}{file}/featurecounts_DE_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno  = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DE", COUNTENV)['OPTIONS'][0].items()),
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
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edgerDE:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
    output: allM= expand("{outdir}/Figures/DE_EDGER_{combi}_DataSet_figure_AllConditionsMDS.png", outdir=outdir, combi=combi),
            allsM= expand("{outdir}/Figures/DE_EDGER_{combi}_DataSet_figure_AllConditionsSumMDS.png", outdir=outdir, combi=combi),
            allBCV =   expand("{outdir}/Figures/DE_EDGER_{combi}_DataSet_figure_AllConditionsBCV.png", outdir=outdir, combi=combi),
            allQLD =   expand("{outdir}/Figures/DE_EDGER_{combi}_DataSet_figure_AllConditionsQLDisp.png", outdir=outdir, combi=combi),
            MDplot =  expand("{outdir}/Figures/DE_EDGER_{combi}_{comparison}_figure_MD.png", outdir=outdir, comparison=compstr, combi=combi),
            allN =   expand("{outdir}/Tables/DE_EDGER_{combi}_DataSet_table_AllConditionsNormalized.tsv.gz", outdir=outdir, combi=combi),
            res =   expand("{outdir}/Tables/DE_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi)
            # dift =  expand("{outdir}/Tables/DE_EDGER_{combi}_{comparison}_EDGER_DE_{comparison}_genes_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["logFC-sorted","pValue-sorted"]),
            # sigdift= expand("{outdir}/Tables/DE_EDGER_{combi}_{comparison}_Sig_EDGER_DE_{comparison}_genes_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["pValue-sorted"]),
    log:    expand("LOGS/{outdir}run_edger.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS,DEBIN]),
            outdir = outdir,
            compare = comparison,
            combi = combi
            ref = ANNOTATION
            # cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DTU')['OPTIONS'][2].items())
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.outdir} {params.combi} {params.compare} {threads} 2> {log} "

rule filter_significant_edgerDE:
    input:  dift = rules.run_edgerDE.output.res
    output: sig = expand("{outdir}/Tables/Sig_DE_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi),
            sig_u = expand("{outdir}/Tables/SigUP_DE_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi),
            sig_d = expand("{outdir}/Tables/SigDOWN_DE_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi)
    log:    expand("LOGS/{outdir}filter_edgerDE.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params: pv_cut = re.findall("\d+\.\d+", get_cutoff_as_string(config, 'DTU').split("-")[0]),
            lfc_cut = re.findall("\d+\.\d+", get_cutoff_as_string(config, 'DTU').split("-")[1])
    shell: "set +o pipefail; for i in {outdir}EDGER_DE*results.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip > {outdir}Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip > {outdir}SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip > {outdir}SigDOWN_$fn;else touch {outdir}Sig_$fn {outdir}SigUP_$fn {outdir}SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_edgerDE.output.allM,
            rules.run_edgerDE.output.allsM,
            rules.run_edgerDE.output.allBCV,
            rules.run_edgerDE.output.allQLD,
            rules.run_edgerDE.output.MDplot,
            rules.run_edgerDE.output.allN,
            rules.run_edgerDE.output.res
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u,
    output: rules.themall.input.Rmd
    log:    expand("LOGS/{outdir}create_summary_snippet.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
