DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string(config,'DAS')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  dendrogram = expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.pdf", combo=combo, scombo=scombo, comparison=compstr),
            csv = expand("DAS/{combo}/Tables/DAS_DIEGO_{scombo}_{comparison}_table_results.csv", combo=combo, scombo=scombo, comparison=compstr),
            sig = expand("DAS/{combo}/Tables/Sig_DAS_DIEGO_{scombo}_{comparison}_table_results.csv", combo=combo, scombo=scombo, comparison=compstr),
            # sig_d = expand("DAS/{combo}/Tables/SigDOWN_DAS_DIEGO_{scombo}_{comparison}_table.csv", combo=combo, scombo=scombo, comparison=compstr),
            # sig_u = expand("DAS/{combo}/Tables/SigUP_DAS_DIEGO_{scombo}_{comparison}_table.csv", combo=combo, scombo=scombo, comparison=compstr),
            Rmd = expand("REPORTS/SUMMARY/RmdSnippets/{combo}.Rmd", combo=combo)


rule featurecount_unique:
    input:  reads = expand("MAPPED/{scombo}/{{file}}_mapped_sorted_unique.bam", scombo=scombo)
    output: tmp   = temp("DAS/{combo}/Featurecounts/{file}_tmp.counts"),
            tmph = temp("DE/{combo}/Featurecounts/{file}_tmp.head.gz"),
            tmpc = temp("DE/{combo}/Featurecounts/{file}_tmp.count.gz"),
            cts   = "DAS/{combo}/Featurecounts/{file}_mapped_sorted_unique.counts.gz"
    log:    "LOGS/DAS/{combo}/{file}_featurecounts_diego_unique.log"
    conda:  ""+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno = ANNOTATION,
            cpara = lambda wildcards: tool_params(wildcards.file, None, config, "DAS", DASENV.split('_')[0])['OPTIONS'].get('COUNT', ""),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else '',
            sortmem = lambda wildcards, threads:  int(30/MAXTHREAD*threads)
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} |gzip > {output.tmph} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S {params.sortmem}% -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> {output.tmpc} && zcat {output.tmph} {output.tmpc} |gzip > {output.cts} && mv {output.tmp}.summary {output.cts}.summary"


rule create_samplemaps:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output: smap = "DAS/{combo}/Tables/samplemap.txt",
            cmap = "DAS/{combo}/Tables/groupings.txt"
    log:    "LOGS/DAS/{combo}/create_samplemaps.log"
    conda:  ""+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd, config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd, config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"


rule prepare_junction_usage_matrix:
    input:  smap = rules.create_samplemaps.output.smap,
            cnd  = expand(rules.featurecount_unique.output.cts, combo=combo, file=samplecond(SAMPLES, config))
    output: tbl = "DAS/{combo}/Tables/{scombo}_junction_table_dexdas.txt.gz",
            anno = "DAS/{combo}/Tables/{scombo}_ANNOTATION.gz"
    log:    "LOGS/DAS/{combo}/prepare_{scombo}_junction_usage_matrix.log"
    conda:  ""+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            dereps = lambda wildcards, input: get_reps(input.cnd, config,'DAS'),
    shell:  "{params.bins}/Analysis/DAS/FeatureCounts2DIEGO.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"


rule create_contrast_files:
    input:  anno = expand(rules.prepare_junction_usage_matrix.output.anno, combo=combo, scombo=scombo)
    output: contrast = expand("DAS/{combo}/Tables/{scombo}_{comparison}_contrast.txt", combo=combo, scombo=scombo, comparison=compstr)
    log:    expand("LOGS/DAS/{combo}/create_contrast_files.log", combo=combo)
    conda:  ""+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            compare=comparison,
            outdir = 'DAS/'+combo+'/Tables',
            pcombo = scombo if scombo != '' else 'none'
    shell:  "python3 {params.bins}/Analysis/DAS/diego_contrast_files.py -a <(zcat {input.anno}) -b {params.pcombo} -c {params.compare} -o {params.outdir} 2> {log}"


rule run_diego:
    input:  tbl = expand(rules.prepare_junction_usage_matrix.output.tbl, combo=combo, scombo=scombo),
            contrast = expand(rules.create_contrast_files.output.contrast, combo=combo, scombo=scombo)
    output: dendrogram = rules.themall.input.dendrogram,
            csv = rules.themall.input.csv
    log:    expand("LOGS/DAS/{combo}_{scombo}_{comparison}/run_diego.log", combo=combo, comparison=compstr, scombo=scombo)
    conda:  ""+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS, DASBIN]),
            dpara = lambda x: tool_params(samplecond(SAMPLES, config)[0], None, config, "DAS", DASENV.split('_')[0])['OPTIONS'].get('DAS', ""),
            compare = compstr,
            outfile = [i.replace(".pdf","") for i in expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.pdf", combo=combo, scombo=scombo, comparison=compstr)]
    shell:  "set +euo pipefail; arr=({input.contrast}); orr=({params.outfile}); orrt=({log}); for i in ${{!arr[@]}}; do basecond=$(head -n 1 ${{arr[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{arr[$i]}} -x $basecond -e -f ${{orr[$i]}} &> ${{orrt[$i]}};done && arr=({input.contrast}); orr=({output.csv}); orrt=({log}); for i in ${{!arr[@]}}; do basecond=$(head -n 1 ${{arr[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{arr[$i]}} -x $basecond {params.dpara} 1> ${{orr[$i]}} 2>> ${{orrt[$i]}};done"


rule filter_significant:
    input:  csv = rules.run_diego.output.csv
    output: sig = rules.themall.input.sig
    log:    "LOGS/DAS/filter_diegoDAS.log"
    conda:  ""+DASENV+".yaml"
    threads: 1
    params: pv_cut = get_cutoff_as_string(config, 'DAS', 'pvalue'),
            lfc_cut = get_cutoff_as_string(config, 'DAS', 'lfc')
    shell: "set +o pipefail; arr=({input.csv}); orr=({output.sig}); for i in \"${{!arr[@]}}\"; do a=\"${{arr[$i]}}\"; fn=\"${{a##*/}}\"; if [[ -s \"$a\" ]];then cat $a|head -n1 > \"${{orr[$i]}}\"; cat $a| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!$F[10]);if ($F[10] eq \"yes\") {{print}}' >> \"${{orr[$i]}}\" &>> {log}; else touch \"${{orr[$i]}}\"; fi; done"


rule convertPDF:
    input: rules.run_diego.output.dendrogram
    output: dendrogram = expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.png", combo=combo, scombo=scombo, comparison=compstr)
    log:    expand("LOGS/DAS/{combo}_{scombo}_{comparison}/convertPDF.log", combo=combo, comparison=compstr, scombo=scombo)
    conda:  ""+DASENV+".yaml"
    threads: MAXTHREAD
    shell: "for pdfile in {input} ; do convert -verbose -density 500 -resize '800' $pdfile ${{pdfile%pdf}}png; done"


rule create_summary_snippet:
    input:  rules.convertPDF.output.dendrogram,
            rules.themall.input.csv,
            rules.themall.input.sig
            #rules.run_diego.output.sig_d,
            #rules.run_diego.output.sig_u
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DAS/{combo}/create_summary_snippet.log", combo=combo)
    conda:  ""+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS,
            abspathfiles = lambda w, input: [os.path.abspath(x) for x in input]
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {params.abspathfiles} --output {output} --env {DASENV} --loglevel DEBUG 2>> {log}"
