rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{src}}/{{dir}}/{{gen}}{{name}}_{map}.idx", ref=REFERENCE, map=MAPPERBIN),
            tmp = temp("TEMP/{src}/{dir}/{gen}{name}.fa"),
            tmpa = temp("TEMP/{src}/{dir}/{gen}{name}.anno")
    log:    expand("LOGS/{{src}}/{{dir}}/{{gen}}{{name}}_{map}.idx.log", map=MAPPERBIN)
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in index_params(str.join(os.sep,[wildcards.dir,wildcards.src]), config, "MAPPING")[0].items()),
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(REFERENCE,wildcards.dir,config["ANNOTATION"][wildcards.dir])),
            genpath = lambda wildcards: "{refe}/{genname}".format(refe=REFERENCE,genname=os.path.join(wildcards.src,wildcards.dir))
    shell: "zcat {input.fa} > {output.tmp} && zcat {params.anno} > {output.tmpa} && rm -rf STARTMP && {params.mapp} --runThreadN {threads} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {output.tmpa} --sjdbGTFtagExonParentTranscript Parent  2> {log} && cat LOG.out >> {log} && rm -f LOG.out"

rule mapping:
    input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
            index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], src = str.join(os.sep, source_from_sample(wildcards.file).split(os.sep)[1:]), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERBIN),
            ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
            tmp = temp("STAROUTTMP/{file}"),
            tmpdir = temp("STARTMP/{file}")
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[1].items()),
            mapp=MAPPERBIN,
            genpath = lambda wildcards: "{ref}/{gen}".format(ref=REFERENCE,gen=genomepath(wildcards.file,config)),
            anno = lambda wildcards: anno_from_file(wildcards.file, config)
    shell: "{params.mapp} {params.mpara}  --runThreadN {threads} --genomeDir {params.genpath} q --readFilesCommand zcat --readFilesIn {input.r1} {input.r2} --outFileNamePrefix {output.tmp} 2>> {log} && cd {output.tmp};cat Aligned.out.sam | tee >(grep -v -P '\t4\t' > {output.mapped}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log}; cd .."
