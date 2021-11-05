MAPPERBIN, MAPPERENV = env_bin_from_config3(config,'MAPPING')

rule generate_index:
    input:  ref = REFERENCE
    output: idx = directory(INDEX),
            uidx = expand("{refd}/INDICES/{mape}_{unikey}/{pref}", refd=REFDIR, mape=MAPPERENV, unikey=get_dict_hash(subDict(tool_params(SAMPLES[0], None, config, "MAPPING", MAPPERENV)["OPTIONS"], ["INDEX"])), pref=PREFIX)
    log:    expand("LOGS/{sets}/{mape}.idx.log", sets=SETS, mape=MAPPERENV)
    conda:  ""+MAPPERENV+".yaml"
    threads: 1
    params: indexer = MAPPERBIN.split(' ')[0],
            ipara = lambda wildcards, input: tool_params(SAMPLES[0], None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('INDEX', ""),
            linkidx = lambda wildcards, output: str(os.path.abspath(str.join(os.sep, str(output.uidx[0]).split(os.sep)[:-1]))) if PREFIX != '' else str(os.path.abspath(str(output.uidx[0]))),
            tolink = lambda wildcards, output: str(os.path.abspath(str.join(os.sep, str(output.idx).split(os.sep)[:-1])))
    shell:  "if [[ -f \"{output.uidx}\/*\" ]]; then ln -fs {params.linkidx} {output.idx} && touch {output.uidx} && echo \"Found bwa index, continue with mapping\" ; else {params.indexer} index -p {output.uidx} {params.ipara} {input.ref} 2> {log} && ln -fs {params.linkidx} {output.idx} && touch {output.uidx};fi"

bwaalg = MAPPERBIN.split(' ')[1]

if bwaalg == 'mem' or MAPPERBIN == 'bwa-mem2':
    if paired == 'paired':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    index = rules.generate_index.output.uidx,
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get("MAP", ""),
                    mapp = MAPPERBIN
                    #idx = lambda wildcards, input: str.join(os.sep,[str(input.index), PREFIX]) if PREFIX != '' else input.index
            shell: "{params.mapp} {params.mpara} -t {threads} {input.index} {input.r1} {input.r2} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

    else:
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    uidx = rules.generate_index.output.uidx[0],
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN
                    #idx = lambda wildcards, input: str.join(os.sep,[str(input.index), PREFIX]) if PREFIX != '' else input.index
            shell:  "{params.mapp} {params.mpara} -t {threads} {input.uidx} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'aln': # not supported as stand alone as we need mappign files to continue the workflow
    if paired == 'paired': # handled like sampe
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    ref = REFERENCE
            output: sai1 = report("MAPPED/{combo}/{file}_mapped.R1.sai", category="MAPPING"),
                    sai2 = report("MAPPED/{combo}/{file}_mapped.R2.sai", category="MAPPING"),
                    mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN,
                    mapp1 = MAPPERBIN.split(' ')[0]
            shell:  "{params.mapp} {params.mpara} {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2}| tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
### FOR LATER IF WE EVER NEED aln MODE
#        rule mapping:
#            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
#                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
#                    #index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file, config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), mape=MAPPERENV, extension=check_tool_params(wildcards.file, None, config, 'MAPPING', 2)),
#                    ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file, config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
#            output: sai1 = report("MAPPED/{combo}/{file}_mapped.R1.sai", category="MAPPING"),
#                    sai2 = report("MAPPED/{combo}/{file}_mapped.R2.sai", category="MAPPING"),
#                    mapped = "UNMAPPED/{combo}/{file}_mapped.sam",
#                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
#            log:    "LOGS/{combo}/{file}/mapping.log"
#            conda:  ""+MAPPERENV+".yaml"
#            threads: MAXTHREAD
#            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'MAPPING')['OPTIONS'][1].items()),
#                    mapp=MAPPERBIN
#            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.r1} > {output.sai1} 2>> {log} && {params.mapp} {params.mpara} -t {threads} {input.ref} {input.r2} > {output.sai2} 2>> {log} && touch {output.unmapped} && touch {output.mapped}"

    else: #handled like samse
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    ref = REFERENCE
            output: sai = report("MAPPED/{combo}/{file}_mapped.sai", category="MAPPING"),
                    mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN,
                    mapp1 = MAPPERBIN.split(' ')[0]
            shell:  "{params.mapp1} aln {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} && {params.mapp} {params.mpara} {input.ref} {output.sai} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
### FOR LATER IF WE EVER NEED aln MODE
#        rule mapping:
#            input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
#                    #index = lambda wildcards: expand(rules.generate_index.output.idx, ref=REFERENCE, dir=source_from_sample(wildcards.file, config), gen=genome(wildcards.file, config), name=namefromfile(wildcards.file, config), mape=MAPPERENV, extension=check_tool_params(wildcards.file, None, config, 'MAPPING', 2)),
#                    ref = lambda wildcards: expand(rules.generate_index.input.fa, ref=REFERENCE, dir = source_from_sample(wildcards.file, config), gen =genome(wildcards.file, config), name=namefromfile(wildcards.file, config))
#            output: sai = report("MAPPED/{combo}/{file}_mapped.sai", category="MAPPING"),
#                    mapped = "UNMAPPED/{combo}/{file}_mapped.sam",
#                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
#            log:    "LOGS/{combo}/{file}/mapping.log"
#            conda:  ""+MAPPERENV+".yaml"
#            threads: MAXTHREAD
#            params: mpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, 'MAPPING')['OPTIONS'][1].items()),
#                    mapp=MAPPERBIN
#            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} 2>> {log} && touch {output.unmapped} && touch {output.mapped}"

elif bwaalg == 'samse':
    rule mapping:
        input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                ref = REFERENCE
        output: sai = report("MAPPED/{combo}/{file}_mapped.sai", category="MAPPING"),
                mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN,
                mapp1 = MAPPERBIN.split(' ')[0]
        shell:  "{params.mapp1} aln {params.mpara} -t {threads} {input.ref} {input.query} > {output.sai} && {params.mapp} {params.mpara} {input.ref} {output.sai} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'sampe':
    rule mapping:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                ref = REFERENCE
        output: sai1 = report("MAPPED/{combo}/{file}_mapped.R1.sai", category="MAPPING"),
                sai2 = report("MAPPED/{combo}/{file}_mapped.R2.sai", category="MAPPING"),
                mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
        log:    "LOGS/{combo}/{file}/mapping.log"
        conda:  ""+MAPPERENV+".yaml"
        threads: MAXTHREAD
        params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                mapp = MAPPERBIN,
                mapp1 = MAPPERBIN.split(' ')[0]
        shell:  "{params.mapp} {params.mpara} {input.ref} {input.sai1} {input.sai2} {input.r1} {input.r2}| tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

elif bwaalg == 'bwasw':
    if paired == 'paired':
        rule mapping:
            input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_R1_trimmed.fastq.gz",
                    r2 = "TRIMMED_FASTQ/{combo}/{file}_R2_trimmed.fastq.gz",
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN
            shell: "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.r1} {input.r2} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"

    else:
        rule mapping:
            input:  query = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz",
                    ref = REFERENCE
            output: mapped = temp(report("MAPPED/{combo}/{file}_mapped.sam", category="MAPPING")),
                    unmapped = "UNMAPPED/{combo}/{file}_unmapped.fastq.gz"
            log:    "LOGS/{combo}/{file}/mapping.log"
            conda:  ""+MAPPERENV+".yaml"
            threads: MAXTHREAD
            params: mpara = lambda wildcards: tool_params(wildcards.file, None, config, 'MAPPING', MAPPERENV)['OPTIONS'].get('MAP', ""),
                    mapp = MAPPERBIN
            shell:  "{params.mapp} {params.mpara} -t {threads} {input.ref} {input.query} | tee >(samtools view -h -F 4 > {output.mapped}) >(samtools view -h -f 4 |samtools fastq -n - | pigz > {output.unmapped}) 1>/dev/null 2>> {log} && touch {output.unmapped}"
