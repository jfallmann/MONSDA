rule trimgalore_trim:
    input:  expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES)
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} --no_report_file --gzip {params.tpara} -o {params.odir} {input} &> {log}"

rule trimgalore_rename:
    input:  rules.trimgalore_trim.output
    output: "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    shell:  "mv {input[0]} {output[0]}"
