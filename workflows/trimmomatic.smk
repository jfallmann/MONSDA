rule trimmomatic_trim:
        input:  "FASTQ/{rawfile}.fastq.gz"
        output: "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz"
        log:    "LOGS/{rawfile}_trimmed.log"
        threads: 1
        conda: "../envs/trimmomatic.yaml"
        shell:  "java -jar {BINS}/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred33 {input[0]} {output} -threads {threads} -trimlog {log} MINLEN:25"
