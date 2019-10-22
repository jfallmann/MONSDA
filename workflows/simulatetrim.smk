if paired == 'paired':
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
        threads: 1
        shell:  "ln -s {input.r1} {output.r1} && ln -s {input.r2} {output.r2}"

else:
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
        output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
        threads: 1
        shell:  "ln -s {input.r1} {output.r1}"
