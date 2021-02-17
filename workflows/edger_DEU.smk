DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

combi = "COMBINATION"

outdir = "DEU/EDGER"
comparison = comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("{outdir}/DEU_EDGER_{combi}_SESSION.gz", outdir=outdir, combi=combi),
            Rmd = "REPORTS/SUMMARY/RmdSnippets/SUM_DEU_EDGER.Rmd"

# rule featurecount_unique:
#     input:  reads = "MAPPED/{combo}{file}_mapped_sorted_unique.bam"
#     output: tmp   = temp(expand("{outdir}/Featurecounts_DEU_edger/{{combo}}{{file}}_tmp.counts", outdir=outdir)),
#             cts   = "DEU/Featurecounts_DEU/{combo}{file}_mapped_sorted_unique.counts"
#     log:    "LOGS/{combo}{file}/featurecounts_DEU_edger_unique.log"
#     conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
#     threads: MAXTHREAD
#     params: countb = COUNTBIN,
#             anno  = ANNOTATION,
#             cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DEU", COUNTENV)['OPTIONS'][0].items()),
#             paired   = lambda x: '-p' if paired == 'paired' else '',
#             stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
#     shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"
#
# rule prepare_count_table:
#     input:   cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
#     output:  tbl  = expand("{outdir}/Tables/COUNTS.gz", outdir=outdir),
#              anno = expand("{outdir}/Tables/ANNOTATION.gz", outdir=outdir)
#     log:     expand("LOGS/{outdir}/prepare_count_table.log", outdir=outdir)
#     conda:   "nextsnakes/envs/"+DEUENV+".yaml"
#     threads: 1
#     params:  dereps = lambda wildcards, input: get_reps(input.cnd, config,'DEU'),
#              bins = BINS,
#     shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edgerDEU:
    input:  tbl  = expand("{outdir}/Tables/{combi}_COUNTS.gz",outdir=outdir, combi=combi),
            anno = expand("{outdir}/Tables/{combi}_ANNOTATION.gz",outdir=outdir, combi=combi)
    # input:  tbl = rules.prepare_count_table.output.tbl,
    #         anno = rules.prepare_count_table.output.anno,
    output: session = rules.themall.input.session,
            allM    = expand("{outdir}/Figures/DEU_EDGER_{combi}_DataSet_figure_AllConditionsMDS.png", outdir=outdir, combi=combi),
            allBCV  = expand("{outdir}/Figures/DEU_EDGER_{combi}_DataSet_figure_AllConditionsBCV.png", outdir=outdir, combi=combi),
            allQLD  = expand("{outdir}/Figures/DEU_EDGER_{combi}_DataSet_figure_AllConditionsQLDisp.png", outdir=outdir, combi=combi),
            MDplot =  expand("{outdir}/Figures/DEU_EDGER_{combi}_{comparison}_figure_MD.png", outdir=outdir, comparison=compstr, combi=combi),
            allN    = expand("{outdir}/Tables/DEU_EDGER_{combi}_DataSet_table_AllConditionsNormalized.tsv.gz", outdir=outdir, combi=combi),
            res =   expand("{outdir}/Tables/DEU_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi)
            # dift = expand("{outdir}/EDGER_DEU_{comparison}_exons_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["logFC-sorted","pValue-sorted"]),
            # sigdift = expand("{outdir}/Sig_EDGER_DEU_{comparison}_exons_{sort}.tsv.gz", outdir=outdir, comparison=compstr, sort=["pValue-sorted"]),
            # plot = expand("{outdir}/EDGER_DEU_{comparison}_MD.png", outdir=outdir, comparison=compstr),
    log:    expand("LOGS/{outdir}/run_edger.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DEUBIN]),
            outdir = outdir,
            combi = combi,
            compare = comparison,
            ref = ANNOTATION
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.outdir} {params.combi} {params.compare} {threads} 2> {log}"

rule filter_significant_edgerDEU:
    input:  dift = rules.run_edgerDEU.output.res
    output: sig     = expand("{outdir}/Tables/Sig_DEU_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi),
            sig_u   = expand("{outdir}/Tables/SigUP_DEU_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi),
            sig_d   = expand("{outdir}/Tables/SigDOWN_DEU_EDGER_{combi}_{comparison}_table_results.tsv.gz", outdir=outdir, comparison = compstr, combi=combi)
    log:    expand("LOGS/{outdir}/filter_edgerDEU.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DEU', 'pval'),
            lfc_cut = get_cutoff_as_string(config, 'DEU', 'lfc')
    shell: "set +o pipefail; for i in {outdir}/EDGER_DEU*results.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip > {outdir}/Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip > {outdir}/SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip > {outdir}/SigDOWN_$fn; else touch {outdir}/Sig_$fn {outdir}/SigUP_$fn {outdir}/SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_edgerDEU.output.allM,
            # rules.run_edgerDEU.output.allsM,
            rules.run_edgerDEU.output.allBCV,
            rules.run_edgerDEU.output.allQLD,
            rules.run_edgerDEU.output.allN,
            rules.run_edgerDEU.output.MDplot,
            rules.run_edgerDEU.output.res
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u,
    output: rules.themall.input.Rmd
    log:    expand("LOGS/{outdir}/create_summary_snippet.log",outdir=outdir)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
