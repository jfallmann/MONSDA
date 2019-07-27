rule bbduk_trim_paired:
    input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
            r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
    output: o1 = "TRIMMED_FASTQ/{file}_r1_val_1.fq.gz",
            o2 = "TRIMMED_FASTQ/{file}_r2_val_2.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: MAXTHREAD
    params: odir=lambda wildcards,output:os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} t={threads} in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={params.ada} {params.tpara}"

rule bbduk_rename_paired:
    input:  o1 = rules.bbduk_trim_paired.output.o1,
            o2 = rules.bbduk_trim_paired.output.o2
    output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"
