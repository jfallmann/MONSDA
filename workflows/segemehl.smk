rule generate_index:
    input: fa = expand("{ref}/{{org}}/{{gen}}.fa.gz", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: idx = "{ref}/{org}/{gen}.idx"
    log:    "LOGS/{ref}/{org}/{gen}_idx.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp=MAPPERBIN
    shell: "{params.mapp} -d {input.fa} -x {output.idx} --threads {threads} 2> {log}"

rule mapping:
    input:  "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1]))
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config)[0].items()),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} -d {params.ref} -i {input.index} -q {input[0]} --threads {threads} -o {output[0]} -u {output[1]} 2> {log}"
