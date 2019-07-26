rule cutadapt_trim:
    input:  "FASTQ/{rawfile}.fastq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
    log:    "LOGS/{rawfile}_trimmed.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/4)
    params: ada=ADAPTERS,
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} -a file:{params.ada} {params.tpara} --cores {threads} -o {output} {input} > {log}"
