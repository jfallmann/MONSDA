rule cutadapt_trim:
    input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(MAXTHREAD/8)
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN,
    shell:  "{params.trim} -a file:{params.ada} {params.tpara} --cores {threads} -o {output.o1} {input.r1} > {log}"

rule cutadapt_rename:
    input:  o1 = rules.cutadapt_trim.output.o1
    output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    shell:  "mv {input.o1} {output.r1}"
