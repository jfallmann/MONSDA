DEBIN, DEENV = env_bin_from_config3(config,'DE')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DE')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("DE/{combo}/DE_DESEQ2_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            pca  = expand("DE/{combo}/Figures/DE_DESEQ2_{scombo}_DataSet_figure_PCA.png", combo=combo, scombo=scombo),
            rld  = expand("DE/{combo}/Tables/DE_DESEQ2_{scombo}_DataSet_table_rld.tsv.gz", combo=combo, scombo=scombo),
            vsd  = expand("DE/{combo}/Tables/DE_DESEQ2_{scombo}_DataSet_table_vsd.tsv.gz", combo=combo, scombo=scombo),
            tbl  = expand("DE/{combo}/Tables/DE_DESEQ2_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            plot = expand("DE/{combo}/Figures/DE_DESEQ2_{scombo}_{comparison}_figure_MA.png", combo=combo, comparison=compstr, scombo=scombo),
            vst  = expand("DE/{combo}/Figures/DE_DESEQ2_{scombo}_DataSet_figure_VST-and-log2.png", combo=combo, scombo=scombo),
            heat = expand("DE/{combo}/Figures/DE_DESEQ2_{scombo}_DataSet_figure_heatmap{i}.png", combo=combo,i=[1,2,3,"-samplebysample"], scombo=scombo),
            heats = expand("DE/{combo}/Figures/DE_DESEQ2_{scombo}_DataSet_figure_heatmap-samplebysample.png", combo=combo,i=[1,2,3,"-samplebysample"], scombo=scombo),
            sig   = expand("DE/{combo}/Tables/Sig_DE_DESEQ2_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_d = expand("DE/{combo}/Tables/SigDOWN_DE_DESEQ2_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            sig_u = expand("DE/{combo}/Tables/SigUP_DE_DESEQ2_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, comparison=compstr, scombo=scombo),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)

rule featurecount_unique:
    input:  reads = expand("MAPPED/{scombo}/{{file}}_mapped_sorted_unique.bam", scombo=scombo)
    output: tmp   = temp("DE/{combo}/Featurecounts/{file}_tmp.counts"),
            cts   = "DE/{combo}/Featurecounts/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/DE/{combo}/{file}_featurecounts_deseq2_unique.log"
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

rule run_deseq2:
    input:  cnt  = expand(rules.prepare_count_table.output.tbl, combo=combo, scombo=scombo),
            anno = expand(rules.prepare_count_table.output.anno, combo=combo, scombo=scombo),
    output: session = rules.themall.input.session,
            pca  = rules.themall.input.pca,
            rld  = rules.themall.input.rld,
            vsd  = rules.themall.input.vsd,
            tbl  = rules.themall.input.tbl,
            plot = rules.themall.input.plot,
            vst  = rules.themall.input.vst,
            heat = rules.themall.input.heat,
            heats = rules.themall.input.heats
    log:    expand("LOGS/DE/{combo}_{scombo}_{comparison}/run_deseq2.log", combo=combo, comparison=compstr, scombo=scombo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = str.join(os.sep, [BINS, DEBIN]),
            outdir = 'DE/'+combo,
            compare = comparison,
            pcombo = scombo if scombo != '' else 'none',
            ref = ANNOTATION,
            depara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DE", DEENV.split('_')[0])['OPTIONS'][1].items()),
    shell:  "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {params.ref} {params.outdir} {params.compare} {params.pcombo} {threads} {params.depara} 2> {log}"

rule filter_significant:
    input:  tbl = rules.run_deseq2.output.tbl
    output: sig = rules.themall.input.sig,
            sig_d = rules.themall.input.sig_d,
            sig_u = rules.themall.input.sig_u
    log:    expand("LOGS/DE/{combo}_{scombo}_{comparison}/filter_deseq2.log", combo=combo, comparison=compstr, scombo=scombo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DE', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DE', 'lfc')
    shell: "set +o pipefail; for i in {input};do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[3] || !$F[6]);if ($F[6] < {params.pv_cut} && ($F[3] <= -{params.lfc_cut} ||$F[3] >= {params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/Tables/Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[3] || !$F[6]);if ($F[6] < {params.pv_cut} && ($F[3] >= {params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/Tables/SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[3] || !$F[6]);if ($F[6] < {params.pv_cut} && ($F[3] <= -{params.lfc_cut}) ){{print}}' |gzip > DE/{combo}/Tables/SigDOWN_$fn; else touch DE/{combo}/Tables/Sig_$fn DE/{combo}/Tables/SigUP_$fn DE/{combo}/Tables/SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_deseq2.output.pca,
            rules.run_deseq2.output.rld,
            rules.run_deseq2.output.vsd,
            rules.run_deseq2.output.tbl,
            rules.run_deseq2.output.plot,
            rules.run_deseq2.output.vst,
            rules.run_deseq2.output.heat,
            rules.run_deseq2.output.heats,
            # rules.filter_significant.output.sig,
            # rules.filter_significant.output.sig_d,
            # rules.filter_significant.output.sig_u
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DE/{combo}/create_summary_snippet.log",combo=combo)
    conda:  "nextsnakes/envs/"+DEENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2>> {log}"
