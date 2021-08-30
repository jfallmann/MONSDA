DEDUPBIN, DEDUPENV = env_bin_from_config3(config, 'DEDUP')

wlparams = ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items()) if tool_params(SAMPLES[0], None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items() else None

type = ['sorted', 'unique'] if not rundedup else ['sorted', 'unique', 'sorted_dedup', 'sorted_unique_dedup']

#wildcard_constraints:
#    rawfile = '|'.join(list(SAMPLES)),
#    read = "R1|R2",
#    outdir = dedupoutdir

#rule dedupthemall:
#    input: expand("{outdir}{file}_{read}_dedup.fastq.gz", outdir=outdir, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("{outdir}{file}_dedup.fastq.gz", outdir=outdir, file=samplecond(SAMPLES, config))

if paired == 'paired':
    if wlparams:
        rule whitelist:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                    r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
            output: wl = "DEDUP_FASTQ/{combo}/{file}_whitelist",
                    td = temp(directory("TMP/UMIWL/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_whitelist.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} whitelist {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --read2-in={input.r2} --stdout={output.wl}"

        rule extract:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                    r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                    wl = rules.whitelist.output.wl
            output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                    o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                    td = temp(directory("TMP/UMIEX/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_extract.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --error-correct-cell --whitelist={input.wl} --stdin={input.r1} --read2-in={input.r2} --stdout={output.o1} --read2-out={output.o2}"
    else:
        rule extract:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                    r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
            output: o1 = "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                    o2 = "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz",
                    td = temp(directory("TMP/UMIEX/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_extract.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --read2-in={input.r2} --stdout={output.o1} --read2-out={output.o2}"

else:
    if wlparams:
        rule whitelist:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
            output: wl = "DEDUP_FASTQ/{combo}/{file}_whitelist",
                    td = temp(directory("TMP/UMIWL/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_whitelist.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} whitelist {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --stdout={output.wl}"

        rule extract:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                    wl = rules.whitelist.output.wl
            output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                    td = temp(directory("TMP/UMIEX/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_extract.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --error-correct-cell --whitelist={input.wl} --stdin={input.r1} --stdout={output.o1}"

    else:
        rule extract:
            input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
            output: o1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz",
                    td = temp(directory("TMP/UMIEX/{combo}/{file}"))
            log:   "LOGS/{combo}/{file}_dedup_extract.log"
            conda: "NextSnakes/envs/"+DEDUPENV+".yaml"
            threads: 1
            params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
                    dedup = DEDUPBIN
            shell:  "mkdir -p {output.td} && {params.dedup} extract {params.dpara} --temp-dir {output.td} --log={log} --stdin={input.r1} --stdout={output.o1}"

rule dedupbam:
        input:  bam = expand("MAPPED/{{combo}}/{{file}}_{type}.bam", type=type)
        output: bam = report("MAPPED/{combo}/{file}_{type}_dedup.bam", category="DEDUP"),
                td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
        log:    "LOGS/{combo}/{file}_{type}/dedupbam.log"
        conda:  "NextSnakes/envs/"+DEDUPENV+".yaml"
        threads: 1
        priority: 0               # This should be done after all mapping is done
        params: dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][2].items()),
                dedup = DEDUPBIN
        shell: "mkdir -p {output.td} && {params.dedup} dedup {params.dpara} --paired --temp-dir {output.td} --stdin={input.bam} --log={log} --stdout={output.bam} 2>> {log} && samtools index {output.bam} 2>> {log}"