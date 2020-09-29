MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  fa = expand("{ref}/{{dir}}/{{gen}}{{name}}.fa.gz", ref=REFERENCE)
    output: idx = expand("{ref}/{{dir}}/{map}/{{extension}}/{{gen}}{{name}}_{{extension}}/{map}.idx", ref=REFERENCE, map=MAPPERENV),
            tmp = temp("TMP/{dir}/{gen}{name}_{extension}.fa"),
    log:    expand("LOGS/{{dir}}/{{gen}}{{name}}_{{extension}}_{map}.idx.log", map=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: MAXTHREAD
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items())
    shell:  "zcat {input.fa} > {output.tmp} && {params.indexer}-build {params.ipara} -p {threads} {output.tmp} {output.idx} 2> {log} && touch {output.idx}"

hs2alg = MAPPERBIN.split(' ')[1] if ' ' in MAPPERBIN else MAPPERBIN

if paired == 'paired':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2))[0].replace('.idx',''),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,

                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else ''
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {input.index} -1 {input.r1} -2 {input.r2} -S {output.mapped} --un-conc-gz {output.unmapped} 2>> {log} && touch {output.unmapped}"

else:
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), map=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
                ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
        output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                stranded = lambda x: '--rna-strandness F' if stranded == 'fr' else '--rna-strandness R' if stranded == 'rf' else ''
        shell: "{params.mapp} {params.mpara} {params.stranded} -p {threads} -x {input.index} -U {input.query} -S {output.mapped} --un-gz {output.unmapped} 2>> {log} && touch {output.unmapped}"
