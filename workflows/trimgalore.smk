rule trimgalore_trim:
    input:  "FASTQ/{rawfile}.fastq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fq.gz"
    log:    "LOGS/{rawfile}_trim.log"
    conda: "../envs/trimgalore.yaml"
    threads: 1
    params: odir=lambda wildcards,output:os.path.dirname(output[0])
    shell:  "trim_galore --no_report_file --gzip -q 25 --length 8 -e 0.15 -o {params.odir} {input[0]} > {log[0]}"

rule trimgalore_rename:
    input:  "TRIMMED_FASTQ/{rawfile}_trimmed.fq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
    shell:  "mv {input} {output}"
