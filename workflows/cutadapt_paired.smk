rule cutadapt_trim_paired:
    input:  r1 = "FASTQ/{rawfile}_r1.fastq.gz",
            r2 = "FASTQ/{rawfile}_r2.fastq.gz"
    output: r1 = "TRIMMED_FASTQ/{rawfile}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{rawfile}_r2_trimmed.fastq.gz"
    log:    "LOGS/{rawfile}_trimmed.log"
    conda: "../envs/"+TRIMMENV+".yaml"
    threads: int(MAXTHREAD/4)
    params: ada=ADAPTERS,
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMMBIN
    shell:  "{params.trim} -a file:{params.ada} {params.tpara} --cores {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log}"
