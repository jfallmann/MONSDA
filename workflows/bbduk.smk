rule bbduk_trim:
    input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: MAXTHREAD
    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN,
    shell:  "{params.trim} t={threads} in={input.r1} out={output.o1} ref={params.ada} {params.tpara}"

rule bbduk_rename:
    input:  o1 = rules.bbduk_trim.output.o1
    output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    shell:  "mv {input.o1} {output.r1}"
