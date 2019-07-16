rule trimgalore_trim_paired:
    input:  r1 = "FASTQ/{rawfile}_r1.fastq.gz",
            r2 = "FASTQ/{rawfile}_r2.fastq.gz"
    output: r1 = "TRIMMED_FASTQ/{rawfile}_r1_val.fq.gz",
            r2 = "TRIMMED_FASTQ/{rawfile}_r2_val.fq.gz"
    log:    "LOGS/{rawfile}_trim.log"
    conda: "../envs/"+TRIMENV+".yaml"
    threads: 1
    params: odir=lambda wildcards,output:os.path.dirname(output.r1[0]),
            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.rawfile, None ,config)[0].items()),
            trim=TRIMBIN
    shell:  "{params.trim} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input} > {log}"

rule trimgalore_rename_paired:
    input:  "TRIMMED_FASTQ/{rawfile}_r1_val.fq.gz",
            "TRIMMED_FASTQ/{rawfile}_r2_val.fq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_r1_trimmed.fastq.gz",
            "TRIMMED_FASTQ/{rawfile}_r2_trimmed.fastq.gz"
    shell:  "mv {input} {output}"
