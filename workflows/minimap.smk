MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{src}}/{{gen}}{{name}}_{{ksize}}_{map}.idx", ref=REFERENCE, map=MAPPERBIN)
    log:    expand("LOGS/{{src}}/{{dir}}/{{gen}}{{name}}_{map}_{{ksize}}.idx.log", map=MAPPERBIN)
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer=MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in index_params(str.join(os.sep,[wildcards.dir,wildcards.src]), config, "MAPPING")["OPTIONS"][0].items()),
    shell: "{params.indexer} {params.ipara} -t {threads} -d {output.idx} {input.fa} 2> {log}"

rule mapping:
    input:  query = lambda wildcards: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz".format(file=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
            index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], src=str.join(os.sep, source_from_sample(wildcards.file).split(os.sep)[1:]), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERBIN, ksize=''.join(tool_params(wildcards.file, None ,config, "MAPPING")[2])[1:]),
            ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir=source_from_sample(wildcards.file).split(os.sep)[0], gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "../envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "MAPPING")["OPTIONS"][1].items()),
            mapp=MAPPERBIN
    shell: "{params.mapp} -t {threads} {params.mpara} {input.index} {input.ref} {input.query} | tee >(grep -v -P '\t4\t' > {output.mapped}) >(grep -P '[^@|\t4\t]' |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log}"
