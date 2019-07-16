rule bbduk_trim:
    input:  r1 = "FASTQ/{file}_r1.fastq.gz",
            r2 = "FASTQ/{file}_r2.fastq.gz",
    output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",
    log:    "LOGS/{file}_trimmed.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: ada={ADAPTERS},
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config)[0].items()),
            trim = TRIMBIN
    shell:  "{params.trim} in1={input.r1[0]} in2={input.r2[0]} out1={output.r1[0]} out2={output.r2[0]} ref={params.ada} {params.tpara}"
