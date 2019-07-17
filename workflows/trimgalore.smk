rule trimgalore_trim:
    input:  "FASTQ/{rawfile}.fastq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fq.gz"
    log:    "LOGS/{rawfile}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: odir=lambda wildcards,output: os.path.dirname(output[0]),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.rawfile, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} --no_report_file --gzip {params.tpara} -o {params.odir} {input} > {log}"

rule trimgalore_rename:
    input:  "TRIMMED_FASTQ/{rawfile}_trimmed.fq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
    shell:  "mv {input[0]} {output[0]}"
