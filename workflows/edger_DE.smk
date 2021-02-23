DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts', 'countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("DE/{combo}/DE_EDGER_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            Rmd = "REPORTS/SUMMARY/RmdSnippets/SUM_DE_EDGER.Rmd"

rule featurecount_unique:
    input:  reads = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("DE/{combo}/Featurecounts_DE_edger/{{file}}_tmp.counts", combo=combo)),
            cts   = "DE/Featurecounts_DE/{combo}/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{combo}/{file}/featurecounts_DE_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno  = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DE", COUNTENV)['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES,config))
    output: tbl  = expand("DE/{combo}/Tables/{scombo}_COUNTS.gz",combo=combo, scombo=scombo),
            anno = expand("DE/{combo}/Tables/{scombo}_ANNOTATION.gz",combo=combo, scombo=scombo)
    log:    expand("LOGS/DE/{combo}/prepare_count_table.log",combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params: dereps = lambda wildcards, input: get_reps(input.cnd,config,'DE'),
            bins = BINS,
    shell:  "{params.bins}/Analysis/build_count_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edgerDE:
    input:  tbl = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno
    output: session = rules.themall.input.session,
            allM= expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsMDS.png", combo=combo, scombo=scombo),
            allsM= expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsSumMDS.png", combo=combo, scombo=scombo),
            allBCV =   expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsBCV.png", combo=combo, scombo=scombo),
            allQLD =   expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsQLDisp.png", combo=combo, scombo=scombo),
            MDplot =  expand("DE/{combo}/Figures/DE_EDGER_{scombo}_{comparison}_figure_MD.png", combo=combo, comparison=compstr, scombo=scombo),
            allN =   expand("DE/{combo}/Tables/DE_EDGER_{scombo}_DataSet_table_AllConditionsNormalized.tsv.gz", combo=combo, scombo=scombo),
            res =   expand("DE/{combo}/Tables/DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo)
            # dift =  expand("DE/{combo}/Tables/DE_EDGER_{scombo}_{comparison}_EDGER_DE_{comparison}_genes_{sort}.tsv.gz", combo=combo, comparison=compstr, sort=["logFC-sorted","pValue-sorted"]),
            # sigdift= expand("DE/{combo}/Tables/DE_EDGER_{scombo}_{comparison}_Sig_EDGER_DE_{comparison}_genes_{sort}.tsv.gz", combo=combo, comparison=compstr, sort=["pValue-sorted"]),
    log:    expand("LOGS/DE/{combo}/run_edger.log",combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DEBIN]),
            combo = combo,
            compare = comparison,
            scombo = scombo,
            ref = ANNOTATION
            # cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'DE')['OPTIONS'][2].items())
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.combo} {params.scombo} {params.compare} {threads} 2> {log} "

rule filter_significant_edgerDE:
    input:  dift = rules.run_edgerDE.output.res
    output: sig = expand("DE/{combo}/Tables/Sig_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_u = expand("DE/{combo}/Tables/SigUP_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_d = expand("DE/{combo}/Tables/SigDOWN_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo)
    log:    expand("LOGS/DE/{combo}/filter_edgerDE.log",combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DE', 'pval'),
            lfc_cut = get_cutoff_as_string(config, 'DE', 'lfc')
    shell: "set +o pipefail; for i in DE/{combo}/EDGER_DE*results.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/SigDOWN_$fn;else touch DE/{combo}/Sig_$fn DE/{combo}/SigUP_$fn DE/{combo}/SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_edgerDE.output.allM,
            rules.run_edgerDE.output.allsM,
            rules.run_edgerDE.output.allBCV,
            rules.run_edgerDE.output.allQLD,
            rules.run_edgerDE.output.MDplot,
            rules.run_edgerDE.output.allN,
            rules.run_edgerDE.output.res
            rules.filter_significant.output.sig,
            rules.filter_significant.output.sig_d,
            rules.filter_significant.output.sig_u
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DE/{combo}/create_summary_snippet.log",combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
