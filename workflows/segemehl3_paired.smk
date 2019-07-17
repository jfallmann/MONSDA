rule generate_index_paired:
    input: fa = expand("{ref}/{{org}}/{{gen}}.fa.gz", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: idx = "{ref}/{org}/{gen}.idx"
    log:    "LOGS/{ref}/{org}/{gen}_idx.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp=MAPPERBIN
    shell: "{params.mapp} -d {input.fa} -x {output.idx} --threads {threads} 2> {log}"

rule mapping_paired:
    input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
            index = lambda wildcards: "{ref}/{gen}{name}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[1]))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None , config, "MAPPING")[0].items()),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} -d {params.ref} -i {input.index} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"
