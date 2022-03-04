DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DAS')
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
            sig   = expand("DAS/{combo}/Tables/Sig_DAS_EDGER_{scombo}_{comparison}_table_resultsDiffSpliceExonTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_d  = expand("DAS/{combo}/Tables/SigDOWN_DAS_EDGER_{scombo}_{comparison}_table_resultsDiffSpliceExonTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            sig_u  = expand("DAS/{combo}/Tables/SigUP_DAS_EDGER_{scombo}_{comparison}_table_resultsDiffSpliceExonTest.tsv.gz", combo=combo, comparison = compstr, scombo=scombo),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

rule featurecount_unique:
    input:  reads = expand("MAPPED/{scombo}/{{file}}_mapped_sorted_unique.bam", scombo=scombo)
    output: tmp   = temp("DAS/{combo}/Featurecounts/{file}_tmp.counts"),
            tmph = temp("DE/{combo}/Featurecounts/{file}_tmp.head.gz"),
            tmpc = temp("DE/{combo}/Featurecounts/{file}_tmp.count.gz"),
            cts   = "DAS/{combo}/Featurecounts/{file}_mapped_sorted_unique.counts.gz"
    log:    "LOGS/DAS/{combo}/{file}_featurecounts_edger_unique.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "DAS", DASENV.split('_')[0])['OPTIONS'].get('COUNT', ""),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} |gzip > {output.tmph} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.tmpc} && zcat {output.tmph} {output.tmpc} |gzip > {output.cts} && mv {output.tmp}.summary {output.cts}.summary"
    

rule prepare_count_table:
    input:   cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output:  tbl  = "DAS/{combo}/Tables/{scombo}_COUNTS.gz",
             anno = "DAS/{combo}/Tables/{scombo}_ANNOTATION.gz"
    log:     "LOGS/DAS/{combo}/{scombo}_prepare_count_table.log"
    conda:   ""+DASENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd, config, 'DAS'),
             bins = BINS
    shell: "{params.bins}/Analysis/build_count_table.py {params.dereps} --ids --table {output.tbl} --anno {output.anno} --loglevel DEBUG 2> {log}"

rule run_edger:
    input:  tbl  = expand(rules.prepare_count_table.output.tbl, combo=combo, scombo=scombo),
            anno = expand(rules.prepare_count_table.output.anno, combo=combo, scombo=scombo),
    output: session = rules.themall.input.session,
            allM    = rules.themall.input.allM,
            allBCV  = rules.themall.input.allBCV,
            allQLD  = rules.themall.input.allQLD,
            resG    = rules.themall.input.resG,
            list    = rules.themall.input.list,
            resS    = rules.themall.input.resS,
            resE    = rules.themall.input.resE
    log:    expand("LOGS/DE/{combo}/run_edger.log", combo=combo)
    conda:  ""+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DASBIN]),
            outdir = 'DAS/'+combo,
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            ref = ANNOTATION
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.tbl} {params.ref} {params.outdir} {params.pcombo} {params.compare} {threads} 2> {log} "

rule filter_significant_edger:
    input:  sort = rules.themall.input.resE
    output: sig= rules.themall.input.sig,
            sig_d= rules.themall.input.sig_d,
            sig_u= rules.themall.input.sig_u,
    log:    "LOGS/DAS/filter_edgerDAS.log"
    conda:  ""+DASENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DAS', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DAS', 'lfc')
    shell: "set +o pipefail; arr=({input.sort}); orr=({output.sig}); orrt=({output.sig_d}); orrr=({output.sig_u}); for i in \"${{!arr[@]}}\"; do a=\"${{arr[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then zcat $a| head -n1 |gzip > \"${{orr[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orrt[$i]}}\"; cp \"${{orr[$i]}}\" \"${{orrr[$i]}}\"; zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] <= -{params.lfc_cut} ||$F[6] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orr[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] >= {params.lfc_cut}) ){{print}}' |gzip >> \"${{orrr[$i]}}\" && zcat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[3]);if ($F[3] < {params.pv_cut} && ($F[6] <= -{params.lfc_cut}) ){{print}}' |gzip >> \"${{orrt[$i]}}\"; else touch \"${{orr[$i]}}\" \"${{orrt[$i]}}\" \"${{orrr[$i]}}\"; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.themall.input.allM,
            rules.themall.input.allBCV,
            rules.themall.input.allQLD,
            rules.themall.input.resG,
            rules.themall.input.resS,
            rules.themall.input.resE,
            rules.themall.input.list,
            rules.themall.input.sig,
            rules.themall.input.sig_d,
            rules.themall.input.sig_u,
            rules.themall.input.session
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DAS/{combo}/create_summary_snippet.log",combo=combo)
    conda:  ""+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DASENV} --loglevel DEBUG 2>> {log}"
