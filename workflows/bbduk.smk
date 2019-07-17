rule bbduk_trim:
    input:  "FASTQ/{file}.fastq.gz"
    output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    log:    "LOGS/{file}_trimmed.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: ada={ADAPTERS},
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim = TRIMBIN
    shell:  "{params.trim} in={input[0]} out={output} ref={params.ada} {params.tpara}"
