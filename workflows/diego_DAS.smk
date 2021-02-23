DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

comparison = comparable_as_string2(config,'DAS')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  Rmd = expand("REPORTS/SUMMARY/RmdSnippets/SUM_DAS_DIEGO.Rmd")

rule featurecount_unique:
    input:  reads = "MAPPED/{scombo}/{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("DAS/{combo}/Featurecounts_DAS_diego/{{scombo}}/{{file}}_tmp.counts", combo=combo)),
            cts   = "DAS/Featurecounts_DAS/{scombo}/{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{scombo}/{file}/featurecounts_DAS_diego_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno  = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DAS", COUNTENV)['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule create_samplemaps:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, scombo=scombo, file=samplecond(SAMPLES, config))
    output: smap = expand("DAS/{combo}/Tables/samplemap.txt", combo=combo),
            cmap = expand("DAS/{combo}/Tables/groupings.txt", combo=combo)
    log:    expand("LOGS/DAS/{combo}/create_samplemaps.log", combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd, config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd, config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"

rule prepare_junction_usage_matrix:
    input:  smap = rules.create_samplemaps.output.smap,
            cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES, config), scombo=scombo)
    output: tbl = expand("DAS/{combo}/Tables/{scombo}_junction_table_dexdas.txt.gz", combo=combo, scombo=scombo),
            anno = expand("DAS/{combo}/Tables/{scombo}_ANNOTATION.gz", combo=combo, scombo=scombo)
    log:    expand("LOGS/DAS/{combo}/prepare_junction_usage_matrix.log", combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            dereps = lambda wildcards, input: get_reps(input.cnd, config,'DAS'),
    shell:  "{params.bins}/Analysis/DAS/FeatureCounts2DIEGO.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule create_contrast_files:
    input:  anno = rules.prepare_junction_usage_matrix.output.anno
    output: contrast = expand("DAS/{combo}/Tables/{scombo}_{comparison}_contrast.txt", combo=combo, scombo=scombo, comparison=compstr)
    log:    expand("LOGS/DAS/{combo}/create_contrast_files.log", combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            compare=comparison,
            combo=combo+'/Tables',
            scombo=scombo
    shell:  "python3 {params.bins}/Analysis/DAS/diego_contrast_files.py -a <(zcat {input.anno}) -b {scombo} -c {params.compare} -o {params.combo} 2> {log}"

rule run_diego:
    input:  tbl = rules.prepare_junction_usage_matrix.output.tbl,
            contrast = rules.create_contrast_files.output.contrast
    output: dendrogram = expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.pdf", combo=combo, scombo=scombo, comparison=compstr),
            csv = expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_table_table.csv", combo=combo, scombo=scombo, comparison=compstr)
    log:    expand("LOGS/DAS/{combo}_{scombo}_{comparison}/run_diego.log", combo=combo, scombo=scombo, comparison=compstr)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS, DASBIN]),
            dpara = lambda x: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(samplecond(SAMPLES, config)[0], None , config, "DAS", "diego")['OPTIONS'][0].items()),
            combo = combo,
            compare = compstr,
            outfile = [i.replace(".pdf","") for i in expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.", combo=combo, scombo=scombo, comparison=compstr)]
    shell:  "array1=({input.contrast}); array2=({params.outfile}); for i in ${{!array1[@]}}; do basecond=$(head -n 1 ${{array1[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{array1[$i]}} -x $basecond -e -f ${{array2[$i]}} 2>> {log};done && array1=({input.contrast}); array2=({output.csv}); for i in ${{!array1[@]}}; do basecond=$(head -n 1 ${{array1[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{array1[$i]}} -x $basecond > ${{array2[$i]}} 2>> {log};done"

rule convertPDF:
    input: rules.run_diego.output.dendrogram
    output: dendrogram = expand("DAS/{combo}/Figures/DAS_DIEGO_{scombo}_{comparison}_figure_dendrogram.png", combo=combo, scombo=scombo, comparison=compstr)
    log:    expand("LOGS/DAS/{combo}/convertPDF.log", combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    shell: "for pdfile in {input} ; do convert -verbose -density 500 -resize '800' $pdfile ${{pdfile%pdf}}png; done"

rule create_summary_snippet:
    input:  rules.convertPDF.output.dendrogram,
            rules.run_diego.output.csv
    output: rules.themall.input.Rmd
    log:    expand("LOGS/DAS/{combo}/create_summary_snippet.log", combo=combo)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: int(MAXTHREAD-1) if int(MAXTHREAD-1) >= 1 else 1
    params: bins = BINS
    shell:  "python3 {params.bins}/Analysis/RmdCreator.py --files {input} --output {output} --loglevel DEBUG 2> {log}"
