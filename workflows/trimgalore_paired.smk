rule trimgalore_trim_paired:
    input:  r1 = expand("FASTQ/{rawfile}_r1.fastq.gz", rawfile=SAMPLES),
            r2 = expand("FASTQ/{rawfile}_r2.fastq.gz", rawfile=SAMPLES)
    output: o1 = "TRIMMED_FASTQ/{file}_r1_val_1.fq.gz",
            o2 = "TRIMMED_FASTQ/{file}_r2_val_2.fq.gz"
    log:    "LOGS/{file}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: int(20/8)
    params: odir=lambda wildcards,output:os.path.dirname(output.o1),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} --cores {threads} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1[0]} {input.r2[0]}&> {log}"

rule trimgalore_rename_paired:
    input:  rules.trimgalore_trim_paired.output
    output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
    shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"
