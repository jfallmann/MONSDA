DEDUPBIN, DEDUPENV = env_bin_from_config2(SAMPLES,config,'DEDUP')

if paired == 'paired':
    rule whitelist:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: whitelist = "FASTQ/{file}_whitelist"
        log:   "LOGS/{file}_dedup_whitelist.log"
        conda: "nextsnakes/envs/"+DEDUPENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards,output:os.path.dirname(output.o1),
                dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][0].items()),
                dedup=DEDUPBIN
        shell:  "{params.dedup} whitelist {params.dpara} --log={log} --stdin={inpur.r1} --read2-in={input.r2} --stdout={output.whitelist}"

    rule extract:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                wl=rules.whitelist.output.whitelist
        output: o1 = "DEDUP_FASTQ/{file}_R1.fastq.gz",
                o2 = "DEDUP_FASTQ/{file}_R2.fastq.gz"
        log:   "LOGS/{file}_dedup_extract.log"
        conda: "nextsnakes/envs/"+DEDUPENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards,output:os.path.dirname(output.o1),
                tpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")['OPTIONS'][0].items()),
                trim=TRIMBIN
        shell:  "{params.dedup} extract {params.dpara} --log={log} --whitelist={input.wl} --stdin={inpur.r1} --read2-in={input.r2} --stdout={output.o1} --read2-out={output.o2}"

else:
    rule whitelist:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: whitelist = "FASTQ/{file}_whitelist"
        log:   "LOGS/{file}_dedup_whitelist.log"
        conda: "nextsnakes/envs/"+DEDUPENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards,output:os.path.dirname(output.o1),
                dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][0].items()),
                dedup=DEDUPBIN
        shell:  "{params.dedup} whitelist {params.dpara} --log={log} --stdin={inpur.r1} --stdout={output.whitelist}"

    rule extract:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                wl=rules.whitelist.output.whitelist
        output: o1 = "DEDUP_FASTQ/{file}_R1.fastq.gz"
        log:   "LOGS/{file}_dedup_extract.log"
        conda: "nextsnakes/envs/"+DEDUPENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards,output:os.path.dirname(output.o1),
                tpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")['OPTIONS'][0].items()),
                trim=TRIMBIN
        shell:  "{params.dedup} extract {params.dpara} --log={log} --whitelist={input.wl} --stdin={inpur.r1} --stdout={output.o1}"


rule dedupbam:
    input:  bam = rules.sam2bam.output.bam
    output: bam = report("MAPPED/{file}_mapped_sorted_dedup.bam", category="DEDUP")
    log:    "LOGS/{file}/dedupbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: MAXTHREAD
    params: fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file))),
            dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][1].items()),
            dedup=DEDUPBIN
    shell: "{params.dedup} dedup {params.dpara} --stdin={input.bam} --log={log} > {output.bam} 2>> {log}"

rule dedupuniqbam:
    input:  bam = rules.sam2bamuniq.output.uniqbam
    output: bam = report("MAPPED/{file}_mapped_sorted_unique_dedup.bam", category="DEDUP")
    log:    "LOGS/{file}/dedupuniqbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: MAXTHREAD
    params: fn = lambda wildcards: "{fn}".format(fn=str(sample_from_path(wildcards.file))),
            dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "DEDUP")['OPTIONS'][1].items()),
            dedup=DEDUPBIN
    shell: "{params.dedup} dedup {params.dpara} --stdin={input.bam} --log={log} > {output.bam} 2>> {log}"
