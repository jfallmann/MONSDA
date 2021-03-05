DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DAS')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("DAS/{combo}/DAS_EDGER_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            allM    = expand("DAS/{combo}/Figures/DAS_EDGER_{scombo}_DataSet_figure_AllConditionsMDS.png", combo=combo, scombo=scombo),
            allBCV  = expand("DAS/{combo}/Figures/DAS_EDGER_{scombo}_DataSet_figure_AllConditionsBCV.png", combo=combo, scombo=scombo),
            allQLD  = expand("DAS/{combo}/Figures/DAS_EDGER_{scombo}_DataSet_figure_AllConditionsQLDisp.png", combo=combo, scombo=scombo),
            resG    = expand("DAS/{combo}/Tables/DAS_EDGER_{scombo}_{comparison}_table_resultsGeneTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            list    = expand("DAS/{combo}/Figures/DAS_EDGER_{scombo}_{comparison}_list_topSpliceSimes.tsv", combo=combo, comparison = compstr, scombo=scombo),
            resS    = expand("DAS/{combo}/Tables/DAS_EDGER_{scombo}_{comparison}_table_resultsSimesTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            resE    = expand("DAS/{combo}/Tables/DAS_EDGER_{scombo}_{comparison}_table_resultsDiffSpliceExonTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

rule featurecount_unique:
    input:  reads = expand("MAPPED/{scombo}/{{file}}_mapped_sorted_unique.bam", scombo=scombo)
    output: tmp   = temp("DAS/{combo}/Featurecounts/{file}_tmp.counts"),
            cts   = "DAS/{combo}/Featurecounts/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/DAS/{combo}/{file}_featurecounts_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DAS", DASENV.split('_')[0])['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output:  tbl  = "DAS/{combo}/Tables/{scombo}_COUNTS.gz",
             anno = "DAS/{combo}/Tables/{scombo}_ANNOTATION.gz"
    log:     "LOGS/DAS/{combo}/{scombo}_prepare_count_table.log"
    conda:   "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd, config, 'DAS'),
             bins = BINS
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edger:
    input:  tbl  = expand(rules.prepare_count_table.output.tbl, combo=combo, scombo=scombo),
            anno = expand(rules.prepare_count_table.output.anno, combo=combo, scombo=scombo),
    output: session = rules.themall.input.session,
            allM    = rules.themall.input.allM,
            allsM   = rules.themall.input.allsM,
            allBCV  = rules.themall.input.allBCV,
            allQLD  = rules.themall.input.allQLD,
            resG    = rules.themall.input.resG,
            list    = rules.themall.input.list,
            resS    = rules.themall.input.resS,
            resE    = rules.themall.input.resE
    log:    expand("LOGS/DAS/{combo}_{scombo}_{comparison}/run_edger.log", combo=combo, comparison = compstr, scombo=scombo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DASBIN]),
            outdir = 'DAS/'+combo,
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            ref = ANNOTATION
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.outdir} {params.pcombo} {params.compare} {threads} 2> {log} "

# rule filter_significant_edger:
#      input:  dift = rules.run_edger.output.dift
#      output: sigdift = rules.themall.input.sigdift
#      log:    expand("LOGS/DAS/{combo}_{scombo}_{comparison}/filter_edgerDAS.log", combo=combo, comparison=compstr, scombo=scombo)
#      conda:  "nextsnakes/envs/"+DASENV+".yaml"
#      threads: 1
#      params: pv_cut = get_cutoff_as_string(config, 'DAS', 'pval'),
#              lfc_cut = get_cutoff_as_string(config, 'DAS', 'lfc')
#      shell: "set +o pipefail; for i in DAS/{combo}/DAS_EDGER*_results_*.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut} ||$F[2] >= {params.lfc_cut}) ){{print}}' |gzip > DAS/{combo}/Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] >= {params.lfc_cut}) ){{print}}' |gzip > DAS/{combo}/SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[2] || !$F[5]);if ($F[5] < {params.pv_cut} && ($F[2] <= -{params.lfc_cut}) ){{print}}' |gzip > DAS/{combo}/SigDOWN_$fn;else touch DAS/{combo}/Sig_$fn DAS/{combo}/SigUP_$fn DAS/{combo}/SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_edger.output.allM,
            rules.run_edger.output.allsM,
            rules.run_edger.output.allBCV,
            rules.run_edger.output.allQLD,
            rules.run_edger.output.resG,
            rules.run_edger.output.resS,
            rules.run_edger.output.resE,
            rules.run_edger.output.list,
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DAS/{combo}/create_summary_snippet.log",combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
