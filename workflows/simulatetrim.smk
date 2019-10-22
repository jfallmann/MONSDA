if paired == 'paired':
    rule simulate_trim:
        input:  r1 = expand("FASTQ/{rawfile}_r1.fastq.gz", rawfile=SAMPLES),
                r2 = expand("FASTQ/{rawfile}_r2.fastq.gz", rawfile=SAMPLES)
        output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w,input: "{r}".format(r=os.path.abspath(input.r1[0])),
                filetolink2 = lambda w,input: "{r}".format(r=os.path.abspath(input.r2[0]))
        shell:  "ln -s {params.filetolink} {output.r1} && ln -s {params.filetolink2} {output.r2}"

else:
    rule simulate_trim:
        input:  r1 = expand("FASTQ/{rawfile}.fastq.gz", rawfile=SAMPLES)
        output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w,input: "{r}".format(r=os.path.abspath(input.r1[0]))
        shell:  "ln -s {params.filetolink} {output.r1}"
