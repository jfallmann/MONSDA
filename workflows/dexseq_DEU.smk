DEUBIN, DEUENV = env_bin_from_config3(config,'DEU')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DEU')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  session = expand("DEU/{combo}/DEU_DEXSEQ_{scombo}_SESSION.gz", combo=combo, scombo=scombo),
            Rmd = "REPORTS/SUMMARY/RmdSnippets/SUM_DEU_DEXSEQ.Rmd"

rule prepare_deu_annotation:
    input:   anno = ANNOTATION
    output:  countgtf = expand("{countanno}", countanno=ANNOTATION.replace('.gtf','_fc_dexseq.gtf')),
             deugtf   = expand("{deuanno}", deuanno=ANNOTATION.replace('.gtf','_dexseq.gtf'))
    log:     "LOGS/featurecount_dexseq_annotation.log"
    conda:   "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params:  bins = BINS,
             countstrand = lambda x: '-s' if stranded == 'fr' or stranded == 'rf' else ''
    shell:  "{params.bins}/Analysis/DEU/prepare_deu_annotation2.py -f {output.countgtf} {params.countstrand} {input.anno} {output.deugtf} 2>> {log}"

rule featurecount_dexseq_unique:
    input:  reads = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            countgtf = expand(rules.prepare_deu_annotation.output.countgtf, ref=REFDIR, countanno=ANNOTATION.replace('.gtf','_fc_dexseq.gtf')),
            deugtf = expand(rules.prepare_deu_annotation.output.deugtf, ref=REFDIR, deuanno=ANNOTATION.replace('.gtf','_dexseq.gtf'))
    output: tmp   = temp(expand("DEU/{combo}/Featurecounts_DEU_dexseq/{{scombo}}/{{file}}_tmp.counts", combo=combo)),
            cts   = "DEU/Featurecounts_DEU_dexseq/{combo}/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{combo}/{file}/featurecounts_dexseq_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb  = COUNTBIN,
            cpara  = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DEU", COUNTENV)['OPTIONS'][0].items()),
            paired = lambda x: '-p' if paired == 'paired' else '',
            bins   = BINS,
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {input.countgtf}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k2,2 -k3,3n -k4,4n -k1,1 -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule prepare_count_table:
    input:   cnd = expand(rules.featurecount_dexseq_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output:  tbl = expand("DEU/{combo}/Tables/COUNTS.gz", combo=combo),
             anno = expand("DEU/{combo}/Tables/ANNOTATION.gz", combo=combo)
    log:     expand("LOGS/DEU/{combo}/prepare_count_table.log", combo=combo)
    conda:   "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    params:  dereps = lambda wildcards, input: get_reps(input.cnd, config,'DEU'),
             bins = BINS,
    shell: "{params.bins}/Analysis/DEU/build_DEU_table.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule run_dexseq:
    # input:  cnt  = expand("DEU/{combo}/Tables/{scombo}_COUNTS.gz",combo=combo, scombo=scombo),
    #         anno = expand("DEU/{combo}/Tables/{scombo}_ANNOTATION.gz",combo=combo, scombo=scombo),
            # flat = expand("{deuanno}",  deuanno=ANNOTATION.replace('.gtf','_dexseq.gtf'))
    input:  cnt  = rules.prepare_count_table.output.tbl,
            anno = rules.prepare_count_table.output.anno,
            flat = rules.prepare_deu_annotation.output.deugtf
    output: html = expand("DEU/{combo}/DEXSeqReport_{scombo}_{comparison}/DEXSeq_{scombo}_{comparison}.html", combo=combo, scombo=scombo, comparison=compstr),
            plot = expand("DEU/{combo}/Figures/DEU_DEXSEQ_{scombo}_{comparison}_figure_DispEsts.png", combo=combo, scombo=scombo, comparison=compstr),
            plotMA  = expand("DEU/{combo}/Figures/DEU_DEXSEQ_{scombo}_{comparison}_figure_plotMA.png", combo=combo, scombo=scombo, comparison=compstr),
            siglist  = expand("DEU/{combo}/Figures/DEU_DEXSEQ_{scombo}_{comparison}_list_sigGroupsFigures.png", combo=combo, scombo=scombo, comparison=compstr),
            tbl  = expand("DEU/{combo}/Tables/DEU_DEXSEQ_{scombo}_{comparison}_table_results.tsv.gz", combo=combo, scombo=scombo, comparison=compstr),
            session = rules.themall.input.session
    log:    expand("LOGS/DEU/{combo}_{scombo}_{comparison}/run_dexseq.log", combo=combo, scombo=scombo, comparison=compstr)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins   = str.join(os.sep,[BINS, DEUBIN]),
            combo = combo,
            compare = comparison,
            scombo = scombo,
            ref = ANNOTATION
    shell: "Rscript --no-environ --no-restore --no-save {params.bins} {input.anno} {input.cnt} {params.ref} {input.flat} {params.combo} {params.scombo} {params.compare} {threads} 2> {log}"

rule filter_significant_dexseq:
    input:  tbl = rules.run_dexseq.output.tbl
    output: sigtbl  = expand("DEU/{combo}/Sig_DEU_DEXSEQ_{comparison}_results.tsv.gz", combo=combo, comparison=compstr)
    log:    expand("LOGS/DEU/{combo}_{comparison}/filter_dexseq.log", combo=combo, comparison=compstr)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DEU', 'pval'),
            lfc_cut = get_cutoff_as_string(config, 'DEU', 'lfc')
    shell: "set +o pipefail; for i in DEU/{combo}/DEU_DEXSEQ_*results.tsv.gz;do fn=\"${{i##*/}}\"; if [[ -s \"$i\" ]];then zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < {params.pv_cut} && ($F[9] <= -{params.lfc_cut} ||$F[9] >= {params.lfc_cut}) ){{print}}' |gzip > DEU/{combo}/Sig_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < {params.pv_cut} && ($F[9] >= {params.lfc_cut}) ){{print}}' |gzip > DEU/{combo}/SigUP_$fn && zcat $i| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[6] || !$F[9]);if ($F[6] < {params.pv_cut} && ($F[9] <= -{params.lfc_cut}) ){{print}}' |gzip > DEU/{combo}/SigDOWN_$fn;else touch DEU/{combo}/Sig_$fn DEU/{combo}/SigUP_$fn DEU/{combo}/SigDOWN_$fn; fi;done 2> {log}"

rule create_summary_snippet:
    input:  rules.run_dexseq.output.plot,
            rules.run_dexseq.output.tbl,
            rules.run_dexseq.output.siglist,
            rules.run_dexseq.output.plotMA,
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DEU/{combo}/create_summary_snippet.log" ,combo=combo)
    conda:  "nextsnakes/envs/"+DEUENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
