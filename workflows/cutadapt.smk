rule cutadapt_trim:
    input:  "FASTQ/{rawfile}.fastq.gz"
    output: "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
    log:    "LOGS/{rawfile}_trimmed.log"
    params: ada=ADAPTERS
    conda: "../envs/cutadapt.yaml"
    threads: 1
    shell:  "cutadapt -a file:{params.ada} -q25 -M 95 -m 8 -e 0.15 -o {output} {input[0]} > {log}"
