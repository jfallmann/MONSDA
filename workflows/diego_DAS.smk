DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads_de']#env_bin_from_config3(config,'COUNTING') ##PINNING subreads package to version 1.6.4 due to changes in 2.0.1 gene_id length cutoff that interfers

outdir = "DAS/DIEGO/"
comparison = comparable_as_string2(config,'DAS')
compstr = [i.split(":")[0] for i in comparison.split(",")]

rule themall:
    input:  dendrogram = expand("{outdir}{comparison}_dendrogram.pdf", outdir=outdir, comparison=compstr),
            csv = expand("{outdir}{comparison}_table.csv", outdir=outdir, comparison=compstr)

rule featurecount_unique:
    input:  reads = "MAPPED/{combo}{file}_mapped_sorted_unique.bam"
    output: tmp   = temp(expand("{outdir}Featurecounts_DAS_diego/{{file}}_tmp.counts", outdir=outdir)),
            cts   = "DAS/Featurecounts_DAS/{combo}{file}_mapped_sorted_unique.counts"
    log:    "LOGS/{combo}{file}/featurecounts_DAS_diego_unique.log"
    conda:  "nextsnakes/envs/"+COUNTENV+".yaml"
    threads: MAXTHREAD
    params: countb = COUNTBIN,
            anno  = ANNOTATION,
            cpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None , config, "DAS", COUNTENV)['OPTIONS'][0].items()),
            paired   = lambda x: '-p' if paired == 'paired' else '',
            stranded = lambda x: '-s 1' if stranded == 'fr' else '-s 2' if stranded == 'rf' else ''
    shell:  "{params.countb} -T {threads} {params.cpara} {params.paired} {params.stranded} -a <(zcat {params.anno}) -o {output.tmp} {input.reads} 2> {log} && head -n2 {output.tmp} > {output.cts} && export LC_ALL=C; tail -n+3 {output.tmp}|sort --parallel={threads} -S 25% -T TMP -k1,1 -k2,2n -k3,3n -u >> {output.cts} && mv {output.tmp}.summary {output.cts}.summary"

rule create_samplemaps:
    input:  cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES, config))
    output: smap = expand("{outdir}Tables/samplemap.txt", outdir=outdir),
            cmap = expand("{outdir}Tables/groupings.txt", outdir=outdir)
    log:    expand("LOGS/{outdir}create_samplemaps.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: slist = lambda wildcards, input: get_diego_samples(input.cnd, config,'DAS'),
            clist = lambda wildcards, input: get_diego_groups(input.cnd, config,'DAS'),
            bins = BINS
    shell:  "echo \'{params.slist}\' 1> {output.smap} 2>> {log} && echo \'{params.clist}\' 1> {output.cmap} 2>> {log}"

rule prepare_junction_usage_matrix:
    input:  smap = rules.create_samplemaps.output.smap,
            cnd  = expand(rules.featurecount_unique.output.cts, file=samplecond(SAMPLES, config))
    output: tbl = expand("{outdir}Tables/junction_table_dexdas.txt.gz", outdir=outdir),
            anno = expand("{outdir}Tables/ANNOTATION.gz", outdir=outdir)
    log:    expand("LOGS/{outdir}prepare_junction_usage_matrix.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            dereps = lambda wildcards, input: get_reps(input.cnd, config,'DAS'),
    shell:  "{params.bins}/Analysis/DAS/FeatureCounts2DIEGO.py {params.dereps} --table {output.tbl} --anno {output.anno} 2> {log}"

rule create_contrast_files:
    input:  anno = rules.prepare_junction_usage_matrix.output.anno
    output: contrast = expand("{outdir}Tables/{comparison}_contrast.txt", outdir=outdir, comparison=compstr)
    log:    expand("LOGS/{outdir}create_contrast_files.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: 1
    params: bins = BINS,
            compare=comparison,
            outdir=outdir+'Tables/'
    shell:  "python3 {params.bins}/Analysis/DAS/diego_contrast_files.py -a <(zcat {input.anno}) -c {params.compare} -o {params.outdir} 2> {log}"

rule run_diego:
    input:  tbl = rules.prepare_junction_usage_matrix.output.tbl,
            contrast = expand(rules.create_contrast_files.output.contrast, outdir=outdir, comparison=compstr)
    output: dendrogram = rules.themall.input.dendrogram,
            csv = rules.themall.input.csv
    log:    expand("LOGS/{outdir}run_diego.log", outdir=outdir)
    conda:  "nextsnakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    params: bins   = str.join(os.sep,[BINS, DASBIN]),
            dpara = lambda x: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(samplecond(SAMPLES, config)[0], None , config, "DAS", DASENV)['OPTIONS'][1].items()),
            outdir = outdir,
            compare = compstr,
            outfile = [i.replace(".pdf","") for i in rules.themall.input.dendrogram]
    shell:  "array1=({input.contrast}); array2=({params.outfile}); for i in ${{!array1[@]}}; do basecond=$(head -n 1 ${{array1[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{array1[$i]}} -x $basecond -e -f ${{array2[$i]}} 2>> {log};done && array1=({input.contrast}); array2=({output.csv}); for i in ${{!array1[@]}}; do basecond=$(head -n 1 ${{array1[$i]}} | awk \'{{print $1}}\'); {params.bins} -a <(zcat {input.tbl}) -b ${{array1[$i]}} -x $basecond > ${{array2[$i]}} 2>> {log};done"
