rule generate_index:
    input: fa = expand("{ref}/{{gen}}{{name}}.fa.gz", ref=config["REFERENCE"])
    output: idx = expand("{ref}/{{gen}}{{name}}_{{ksize}}.idx", ref=config["REFERENCE"])
    log:    expand("LOGS/{ref}/{{gen}}{{name}}_{{ksize}}_idx.log",ref=config["REFERENCE"])
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp=MAPPERBIN,
            ksize = lambda wildcards: str(wildcards.ksize[1:])
    shell: "{params.mapp} -t {threads} -k {params.ksize} -d {output.idx} {input.fa}  2> {log}"

rule mapping:
    input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            idx = lambda wildcards: expand(rules.generate_index.output.idx,ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config), ksize=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1]))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[0].items()),
            index = lambda wildcards: "{ref}/{gen}{name}_{kmer}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config), kmer=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1])),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN
    shell: "{params.mapp} -t {threads} {params.mpara} {params.index} {params.ref} {input.query} | tee >(grep -v -P '\t4\t' > {output.mapped}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log}"
