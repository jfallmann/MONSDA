DASBIN, DASENV = env_bin_from_config3(config,'DAS')
COUNTBIN, COUNTENV = ['featureCounts','countreads']#env_bin_from_config2(SAMPLES,config,'COUNTING')

outdir="DAS/DIEGO/"
comparison=comparable_as_string2(config,'DAS')
comps = comparison.split(",")
files=samplecond(SAMPLES,config)


rule themall:
    input: tbl = expand("{outdir}DEXSeq_{comparison}.tsv.gz", outdir=outdir, comparison=comparison.split(",")),
           plot = expand("{outdir}DEXSeq_{comparison}_DispEsts.pdf", outdir=outdir, comparison=comparison.split(",")),
           html = expand("{outdir}DEXSeqReport_{comparison}/DEXSeq_{comparison}.html", outdir=outdir, comparison=comparison.split(",")),
           session = expand("{outdir}DEXSeq_SESSION.gz", outdir=outdir)# R object?

rule create_genome_annotation_file:
    input:  lambda wildcards: str.join(os.sep,[config["REFERENCE"],os.path.dirname(genomepath(wildcards.file, config)),tool_params(wildcards.file, None, config, 'DAS')['ANNOTATION']])
    output: expand("{outdir}annotation_DIEGO.bed", outdir=outdir)
    log:    "LOGS/DAS/DIEGO/create_genome_annotation_file.log"
    conda:  "snakes/envs/"+DASENV+".yaml"
    threads: MAXTHREAD
    shell:  "perl scripts/Analysis/DAS/DIEGO/gfftoDIEGObed.pl -g  <(awk 'OFS="\t" {if (NR > 5) $1="chr"$1; print}' <(zcat {input})) -o {output} 2> {log}"

for contrast in comps:

    rule



    rule create_list_of_starjunction_files:
        input:  expand("MAPPED/{files}SJ.out.tab",files=samplecond(SAMPLES,config))
        output: expand("{outdir}names_and_data_list.txt")
        log:    "LOGS/DAS/DIEGO/create_list_of_starjunction_files.log"
        conda:  "snakes/envs/"+DASENV+".yaml"
        threads: MAXTHREAD
        shell:  'for i in {input}; do echo -e "$(echo $i | awk -F/ '{print$NF}' | tr -d SJ.out.tab)\t$i >> {output}"; done'

    if "/star/" in samplecond(SAMPLES,config):

        rule prepare_junction_usage_matrix:
            input:  anno=rule.create_genome_annotation_file.output,
                    list=rule.create_list_of_starjunction_files.output
            output: tbl="junction_table.txt"
            log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
            conda:  "snakes/envs/"+DASENV+".yaml"
            threads: MAXTHREAD
            params: out=outdir
            shell:  "python3 scripts/Analysis/DAS/DIEGO/pre_STAR.py -l {input.list} -d {input.anno} -o {params.out}"

    elif "/segemehl/" in samplecond(SAMPLES,config):

        rule prepare_junction_usage_matrix:
            input:  anno=rule.create_genome_annotation_file.output,
                    list=rule.create_list_of_starjunction_files.output
            output: tbl="junction_table.txt"
            log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
            conda:  "snakes/envs/"+DASENV+".yaml"
            threads: MAXTHREAD
            params: out=outdir
            shell:  "perl scripts/Analysis/DAS/DIEGO/pre_segemehl.py -l {input.list} -d {input.anno} -o {params.out}"

    elif "/htseq/" in samplecond(SAMPLES,config):

        rule prepare_junction_usage_matrix:
            input:  anno=rule.create_genome_annotation_file.output,
                    list=rule.create_list_of_starjunction_files.output
            output: tbl="junction_table.txt"
            log:    "LOGS/DAS/DIEGO/prepare_junction_usage_matrix"
            conda:  "snakes/envs/"+DASENV+".yaml"
            threads: MAXTHREAD
            params: out=outdir
            shell:  "perl scripts/Analysis/DAS/DIEGO/HTseq2DIEGO.pl -l {input.list} -d {input.anno} -o {params.out}"

    rule run_diego:
        input:  tbl=
                anno= rule.prepare_junction_usage_matrix.output,
                list=
        output: expand("{outdir}dendrogram", outdir=outdir)
        log:    "LOGS/"
        conda:  "snakes/envs/"+DASENV+".yaml"
        threads: MAXTHREAD
        params: bins   = str.join(os.sep,[BINS,DASBIN]),
                outdir = outdir,
                compare = comparison
        shell:  "python {params.bins} -a {input.tbl} -b {input.anno} -x your_base_condition -e [-f {output}]"


onsuccess:
    print("Workflow finished, no error")
onerror:
	print("ERROR: "+str({log}))
