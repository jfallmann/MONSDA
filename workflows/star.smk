rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{src}}/{{gen}}{{name}}_{map}.idx", ref=REFERENCE, map=MAPPERBIN),
            tmp = temp("TMPIDX/{src}/{dir}/{gen}{name}.fa"),
            tmpa = temp("TMPIDX/{src}/{dir}/{gen}{name}.anno")
    log:    expand("LOGS/{{src}}/{{dir}}/{{gen}}{{name}}_{map}.idx.log", map=MAPPERBIN)
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in index_params(str.join(os.sep,[wildcards.dir,wildcards.src]), config, "MAPPING")[0].items()),
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(REFERENCE,wildcards.dir,config["ANNOTATION"][wildcards.dir])),
            genpath = lambda wildcards: "{ref}/{dir}/{map}/{src}".format(ref=REFERENCE, dir=wildcards.dir, map=MAPPERBIN, src = wildcards.src)
    shell:  "if [[ -f \"{params.genpath}/SAindex\" ]]; then touch {output.idx} {output.tmp} {output.tmpa} && echo \"Found SAindex, continue with mapping\" ; else zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && rm -rf TMPSTAR && {params.mapp} {params.ipara} --runThreadN {threads} --runMode genomeGenerate --outTmpDir TMPSTAR --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} --sjdbGTFtagExonParentTranscript Parent  2> {log} && touch {output.idx} && cat Log.out >> {log} && rm -f Log.out && rm -rf TMPSTAR ;fi"

rule mapping:
    input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], src=str.join(os.sep, source_from_sample(wildcards.file).split(os.sep)[1:]), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERBIN),
            ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
            tmp = temp("TMPSTAROUT/{file}")
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[1].items()),
            mapp=MAPPERBIN,
            genpath = lambda wildcards: "{ref}/{dir}/{map}/{src}".format(ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], src=str.join(os.sep, source_from_sample(wildcards.file).split(os.sep)[1:]), map=MAPPERBIN),
            anno = lambda wildcards: anno_from_file(wildcards.file, config),
            tocopy = lambda wildcards, output: os.path.dirname(output.mapped)
    shell: "{params.mapp} {params.mpara} --runThreadN {threads} --genomeDir {params.genpath} q --readFilesCommand zcat --readFilesIn {input.query} --outFileNamePrefix {output.tmp} 2>> {log} && cat {output.tmp}Aligned.out.sam | tee >(samtools view -h -F 4 - > {output.mapped}) >(samtools view -h -f 4 - |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && rm {output.tmp}Aligned.out.sam && mv {output.tmp}*.out* {params.tocopy} 2>>{log} && touch {output.tmp} 2>>{log}"
