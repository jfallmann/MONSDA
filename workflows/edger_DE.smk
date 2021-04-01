DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts', 'countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("DE/{combo}/DE_EDGER_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            allM    = expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsMDS.png", combo=combo, scombo=scombo),
            allBCV  = expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsBCV.png", combo=combo, scombo=scombo),
            allQLD  = expand("DE/{combo}/Figures/DE_EDGER_{scombo}_DataSet_figure_AllConditionsQLDisp.png", combo=combo, scombo=scombo),
            MDplot  = expand("DE/{combo}/Figures/DE_EDGER_{scombo}_{comparison}_figure_MD.png", combo=combo, comparison=compstr, scombo=scombo),
            allN    = expand("DE/{combo}/Tables/DE_EDGER_{scombo}_DataSet_table_AllConditionsNormalized.tsv.gz", combo=combo, scombo=scombo),
            res     = expand("DE/{combo}/Tables/DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sort    = expand("DE/{combo}/Tables/DE_EDGER_{scombo}_{comparison}_table_results{sort}.tsv.gz", combo=combo, comparison=compstr, scombo=scombo, sort=["LogFCsorted","PValueSorted"]),
            sig     = expand("DE/{combo}/Tables/Sig_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_u   = expand("DE/{combo}/Tables/SigUP_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_d   = expand("DE/{combo}/Tables/SigDOWN_DE_EDGER_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            Rmd     = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

rule featurecount_unique:
    input:  reads = expand("MAPPED/{scombo}/{{file}}_mapped_sorted_unique.bam", scombo=scombo)
    output: tmp   = temp("DE/{combo}/Featurecounts/{file}_tmp.counts"),
            cts   = "DE/{combo}/Featurecounts/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/DE/{combo}/{file}_featurecounts_edger_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DE", DEENV.split('_')[0])['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output:  tbl  = "DE/{combo}/Tables/{scombo}_COUNTS.gz",
             anno = "DE/{combo}/Tables/{scombo}_ANNOTATION.gz"
    log:     "LOGS/DE/{combo}/{scombo}_prepare_count_table.log"
    conda:   "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd, config, 'DE'),
             bins = BINS
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edger:
    input:  tbl = expand(rules.prepare_count_table.output.tbl, combo=combo, scombo=scombo),
            anno = expand(rules.prepare_count_table.output.anno, combo=combo, scombo=scombo)
    output: session = rules.themall.input.session,
            allM    = rules.themall.input.allM,
            allBCV  = rules.themall.input.allBCV,
            allQLD  = rules.themall.input.allQLD,
            MDplot  = rules.themall.input.MDplot,
            allN    = rules.themall.input.allN,
            res     = rules.themall.input.res,
            sort    = rules.themall.input.sort
    log:    expand("LOGS/DE/{combo}/run_edger.log", combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DEBIN]),
            outdir = 'DE/'+combo,
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            ref = ANNOTATION,
            depara = lambda wildcards: list(tool_params(samplecond(SAMPLES, config)[0], None , config, "DE", DEENV.split('_')[0])['OPTIONS'][1].values())[0] if len(tool_params(samplecond(SAMPLES, config)[0], None , config, "DE", DEENV.split('_')[0])['OPTIONS']) > 1 else ''
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.outdir} {params.compare} {params.pcombo} {threads} {params.depara} 2> {log}"

rule filter_significant:
    input:  sort = rules.run_edger.output.sort
    output: sig = rules.themall.input.sig,
            sig_d = rules.themall.input.sig_d,
            sig_u = rules.themall.input.sig_u
    log:    "LOGS/DE/filter_edgerDE.log"
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DE', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DE', 'lfc')
    shell:  "set +o pipefail; array1=({input.sort}); array2=({output.sig}); array3=({output.sig_d}); array4=({output.sig_u}); for i in ${{!array1[@]}}; do a=${{array1[$i]}}; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > ${{array2[$i]}}; cp ${{array2[$i]}} ${{array3[$i]}}; cp ${{array2[$i]}} ${{array4[$i]}}; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] <= -{params.lfc_cut} ||$F[6] >= {params.lfc_cut}) ){{print}}' |gzip >> ${{array2[$i]}} && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] >= {params.lfc_cut}) ){{print}}' |gzip >> ${{array4[$i]}} && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] <= -{params.lfc_cut}) ){{print}}' |gzip >> ${{array3[$i]}}; else touch ${{array2[$i]}} ${{array3[$i]}} ${{array4[$i]}}; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_edger.output.allM,
            rules.run_edger.output.allBCV,
            rules.run_edger.output.allQLD,
            rules.run_edger.output.MDplot,
            rules.run_edger.output.allN,
            rules.run_edger.output.res,
            rules.run_edger.output.sort#,
            #rules.filter_significant.output.sig,
            #rules.filter_significant.output.sig_d,
            #rules.filter_significant.output.sig_u
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DE/{combo}/create_summary_snippet.log", combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
