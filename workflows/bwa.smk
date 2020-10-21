MAPPERBIN, MAPPERENV = env_bin_from_config2(SAMPLES,config,'MAPPING')

rule generate_index:
    input:  ref = REFERENCE
    output: idx = INDEX,
            uidx = expand("{refd}/INDICES/{mape}/{unikey}/idx", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0]))
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
    threads: 1
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None, config, 'MAPPING')['OPTIONS'][0].items()),
            linkidx = lambda wildcards, output: str(os.path.abspath(str.join(os.sep,[str(output.uidx[0]).split(os.sep)][:-1]))),
    shell:  "{params.indexer} index -p {output.uidx} {params.ipara} {input.ref} 2> {log} && ln -s {params.linkidx} {output.idx}"

bwaalg = MAPPERBIN.split(' ')[1]

if bwaalg == 'mem':
    if paired == 'paired':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                    index = rules.generate_index.output.idx,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN
            shell: "{params.mapp} {params.mpara} -t {threads} {input.index} {input.r1} {input.r2} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

    else:
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                    index = rules.generate_index.output.idx,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN
            shell:  "{params.mapp} {params.mpara} -t {threads} {input.index} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'aln': # not supported as stand alone as we need mappign files to continue the workflow
    if paired == 'paired': # handled like sampe
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                    ref = REFERENCE
            output: sai1 = report("MAPPED/{file}_mapped.R1.sai", category="MAPPING"),
                    sai2 = report("MAPPED/{file}_mapped.R2.sai", category="MAPPING"),
                    mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN,
                    mapp1=MAPPERBIN.split(' ')[0]
            shell:  "{params.mapp} {params.mpara} {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2}| tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
### FOR LATER IF WE EVER NEED aln MODE
#        rule mapping:
#            input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
#                    r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
#                    #index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), mape=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
#                    ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
#            output: sai1 = report("MAPPED/{file}_mapped.R1.sai", category="MAPPING"),
#                    sai2 = report("MAPPED/{file}_mapped.R2.sai", category="MAPPING"),
#                    mapped = "UNMAPPED/{file}_mapped.sam",
#                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
#            log:    "LOGS/{file}/mapping.log"
#            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
#            threads: MAXTHREAD
#            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
#                    mapp=MAPPERBIN
#            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.r1} > {output.sai1} 2>> {log} && {params.mapp} {params.mpara} -t {threads} {input.ref} {input.r2} > {output.sai2} 2>> {log} && touch {output.unmapped} && touch {output.mapped}"

    else: #handled like samse
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                    ref = REFERENCE
            output: sai = report("MAPPED/{file}_mapped.sai", category="MAPPING"),
                    mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN,
                    mapp1=MAPPERBIN.split(' ')[0]
            shell:  "{params.mapp1} aln {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} && {params.mapp} {params.mpara} {input.ref} {output.sai} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
### FOR LATER IF WE EVER NEED aln MODE
#        rule mapping:
#            input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
#                    #index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file,config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), mape=MAPPERENV, extension=check_tool_params(wildcards.file, None ,config, 'MAPPING', 2)),
#                    ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file,config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
#            output: sai = report("MAPPED/{file}_mapped.sai", category="MAPPING"),
#                    mapped = "UNMAPPED/{file}_mapped.sam",
#                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
#            log:    "LOGS/{file}/mapping.log"
#            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
#            threads: MAXTHREAD
#            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
#                    mapp=MAPPERBIN
#            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} 2>> {log} && touch {output.unmapped} && touch {output.mapped}"

elif bwaalg == 'samse':
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                ref = REFERENCE
        output: sai = report("MAPPED/{file}_mapped.sai", category="MAPPING"),
                mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                mapp1=MAPPERBIN.split(' ')[0]
        shell:  "{params.mapp1} aln {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} && {params.mapp} {params.mpara} {input.ref} {output.sai} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'sampe':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                ref = REFERENCE
        output: sai1 = report("MAPPED/{file}_mapped.R1.sai", category="MAPPING"),
                sai2 = report("MAPPED/{file}_mapped.R2.sai", category="MAPPING"),
                mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
        log:    "LOGS/{file}/mapping.log"
        conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                mapp=MAPPERBIN,
                mapp1=MAPPERBIN.split(' ')[0]
        shell:  "{params.mapp} {params.mpara} {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2}| tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'bwasw':
    if paired == 'paired':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz",
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN
            shell: "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.r1} {input.r2} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

    else:
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{file}_unmapped.fastq.gz"
            log:    "LOGS/{file}/mapping.log"
            conda:  "nextsnakes/envs/"+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'MAPPING')['OPTIONS'][1].items()),
                    mapp=MAPPERBIN
            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
