MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx1 = expand("{ref}/{{dir}}/{map}/{{gen}}{{name}}_{{extension}}_{map}.idx", ref=REFERENCE, map=MAPPERENV),
            idx2 = expand("{ref}/{{dir}}/{map}/{{gen}}{{name}}_{{extension}}_{map}_second.idx", ref=REFERENCE, map=MAPPERENV)
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{extension}}_{map}.idx.log", map=MAPPERENV)
    conda:  "snakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN,
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items())
    shell: "{params.indexer} --threads {threads} {params.ipara} -d {input.fa} -x {output.idx} -y {output.idx2} 2> {log}"

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                index2 = lambda wildcards: expand(rules.generate_index.output.idx2, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -j {input.index2} -q {input.r1} -p {input.r2} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                index2 = lambda wildcards: expand(rules.generate_index.output.idx2, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = report("MAPPED/{file}_mapped.sam", category="MAPPING"),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "snakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN
        shell: "{params.mapp} {params.mpara} -d {input.ref} -i {input.index} -j {input.index2} -q {input.query} --threads {threads} -o {output.mapped} -u {output.unmapped} 2> {log}"
