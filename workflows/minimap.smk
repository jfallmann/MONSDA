rule generate_index:
    input: fa = expand("{ref}/{{org}}/{{gen}}.fa.gz", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: idx = expand("{ref}/{org}/{gen}_{ksize}.idx",ksize=lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1])),)
    log:    "LOGS/{ref}/{org}/{gen}_idx.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp=MAPPERBIN,
            ksize = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1])),
    shell: "{params.mapp} -k {params.ksize} -d {input.fa} -x {output.idx} --threads {threads} 2> {log}"

rule mapping:
    input:  expand("TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",rawfile=SAMPLES)
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[0].items()),
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1])),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN
    shell: "{params.mapp} -t {threads} {params.mpara} {params.index} {params.ref} {input[0]} | tee >(grep -v -P '\t4\t' > {output[0]}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output[1]}) 1>/dev/null 2>> {log}"
