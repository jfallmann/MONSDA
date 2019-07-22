rule generate_index:
    input: fa = expand("{ref}/{{gen}}{{name}}.fa.gz", ref=config["REFERENCE"])
    output: idx = expand("{ref}/{{gen}}{{name}}_{map}.idx", ref=config["REFERENCE"], map=MAPPERBIN),
            tmp = temp("TEMP/{gen}{name}.fa")
    log:    expand("LOGS/{ref}/{{gen}}{{name}}_idx.log",ref=config["REFERENCE"])
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp = MAPPERBIN,
            mpara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in idx_params(input.fa, None, config, "MAPPING")[0].items()),
            anno = lambda wildcards: "{annotation}".format(annotation=os.path.join(config["REFERENCE"],os.path.dirname(wildcards.gen),config["ANNOTATION"][os.path.dirname(wildcards.gen)])),
            genpath = lambda wildcards: "{refe}/{genname}".format(refe=config["REFERENCE"],genname=os.path.dirname(wildcards.gen))
    shell: "zcat {input.fa} > {output.tmp} && rm -rf STARTMP && {params.mapp} --runThreadN {threads} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir {params.genpath} --genomeFastaFiles {output.tmp} --sjdbGTFfile {params.anno} --sjdbGTFtagExonParentTranscript Parent  2> {log} && cat LOG.out >> {log} && rm -f LOG.out"

rule mapping:
    input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            index = lambda wildcards: "{ref}/{gen}{name}_{map}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[2]), map=MAPPERBIN)
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
            tmp = temp("STAROUTTMP/{file}")
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[1].items()),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN,
            genpath = lambda wildcards: "{ref}/{gen}".format(ref=REFERENCE,gen=genomepath(wildcards.file,config))
    shell: "{params.mapp} {params.mpara}  --runThreadN {threads} --genomeDir {params.genpath} q --readFilesCommand zcat --readFilesIn {input.query} --outFileNamePrefix {output.tmp} 2> {log} && cd {output.tmp};cat Aligned.out.sam | tee >(grep -v -P '\t4\t' > {output.mapped}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} ;cd .."
