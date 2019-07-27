rule trimgalore_trim:
    input:  r1 = "FASTQ/{rawfile}.fastq.gz"
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/8)
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN,
    shell:  "{params.trim} --cores {threads} --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1} &> {log}"

rule trimgalore_rename:
    input:  o1 = expand(rules.trimgalore_trim.output, rawfile=SAMPLES)
    output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/8)
    shell:  "mv {input.o1} {output.r1}"
