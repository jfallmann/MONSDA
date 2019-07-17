rule trimgalore_trim_paired:
    input:  r1 = expand("FASTQ/{rawfile}_r1.fastq.gz", rawfile=SAMPLES),
            r2 = expand("FASTQ/{rawfile}_r2.fastq.gz", rawfile=SAMPLES)
    output: r1 = "TRIMMED_FASTQ/{file}_r1_val.fq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_val.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: odir=lambda wildcards,output:os.path.dirname(output.r1[0]),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input} > {log}"

#rule trimgalore_trim_paired:
#    input:  r1 = expand("FASTQ/{rawfile}_r1.fastq.gz", rawfile=SAMPLES),
#            r2 = expand("FASTQ/{rawfile}_r2.fastq.gz", rawfile=SAMPLES)
#    output: r1 = expand("TRIMMED_FASTQ/{file}_r1_val.fq.gz",file=samplecond(SAMPLES,config)),
#            r2 = expand("TRIMMED_FASTQ/{file}_r2_val.fq.gz",file=samplecond(SAMPLES,config))
#    log:    expand("LOGS/{file}_trim.log",file=samplecond(SAMPLES,config))
#    conda: "../envs/"+TRIMENV+".yaml"
#    threads: 1
#    params: odir=lambda wildcards,output:os.path.dirname(output.r1[0]),
#            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
#            trim=TRIMBIN
#    shell:  "{params.trim} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input} > {log}"

rule trimgalore_rename_paired:
    input:  rules.trimgalore_trim_paired.output
    output: "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
    shell:  "mv {input} {output}"
