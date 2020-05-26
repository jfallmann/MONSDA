MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')
rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{gen}}{{name}}_{{ksize}}_{map}.idx", ref=REFERENCE, map=MAPPERENV)
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{ksize}}_{map}.idx.log", map=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer=MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
    shell: "{params.indexer} {params.ipara} -t {threads} -d {output.idx} {input.fa} 2> {log}"

rule mapping:
    input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, ksize=tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][2]),
            ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
    output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
            unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
    log:    "LOGS/{file}/mapping.log"
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params:  mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
            mapp=MAPPERBIN
    shell: "{params.mapp}  {params.mpara} -t {threads} {input.index} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
