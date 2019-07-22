rule generate_index_paired:
    input: fa = expand("{ref}/{{org}}/{{gen}}.fa.gz", ref=config["REFERENCE"])#, gen=pathstogenomes(SAMPLES, config))
    output: idx = expand("{{ref}}/{{org}}/{{gen}}_{map}.idx", map=MAPPERBIN)
    log:    "LOGS/{ref}/{org}/{gen}_idx.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp=MAPPERBIN,
            mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[0].items())
    shell: "{params.mapp} {params.mpara} -d {input.fa} -x {output.idx} --threads {threads} 2> {log}"

rule mapping_paired:
    input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
            index = lambda wildcards: "{ref}/{gen}{name}_{map}.idx".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[2]), map=MAPPERBIN)
    output: report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None , config, "MAPPING")[1].items()),
            ref = lambda wildcards: check_ref("{ref}/{gen}{name}.fa".format(ref=REFERENCE,gen=genomepath(wildcards.file,config), name=namefromfile(wildcards.file, config))),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} -d {params.ref} -i {input.index} -q {input.r1[0]} -p {input.r2[0]} --threads {threads} -o {output[0]} -u {output[1]} 2> {log}"
