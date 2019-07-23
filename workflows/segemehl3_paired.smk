rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{src}}/{{gen}}{{name}}_{map}.idx", ref=REFERENCE, map=MAPPERBIN)
    log:    expand("LOGS/{{src}}/{{dir}}/{{gen}}{{name}}_{map}.idx.log", map=MAPPERBIN)
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mapp = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in index_params(str.join(os.sep,[wildcards.dir,wildcards.src]), config, "MAPPING")[0].items())
    shell: "{params.mapp} {params.ipara} -d {input.fa} -x {output.idx} --threads {threads} 2> {log}"

rule mapping_paired:
    input:  r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
            index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], src=str.join(os.sep, source_from_sample(wildcards.file).split(os.sep)[1:]), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERBIN),
            ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file).split(os.sep)[0], gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz",
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: 20
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")[1].items()),
            mapp=MAPPERBIN
    shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"
