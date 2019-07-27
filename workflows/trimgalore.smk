rule trimgalore_trim:
    input:  r1 = expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES)
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/8)
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN,
    shell:  "{params.trim} --cores {threads} --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1} &> {log}"

rule trimgalore_rename:
    input:  r1 = rules.trimgalore_trim.output
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/8)
    shell:  "mv {input.r1} {output.o1}"
