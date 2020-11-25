
if paired == 'paired':
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not rundedup else "DEDUP_FASTQ/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not rundedup else "DEDUP_FASTQ/{file}_R2_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w,input: "{r}".format(r=os.path.abspath(input.r1[0])),
                filetolink2 = lambda w,input: "{r}".format(r=os.path.abspath(input.r2[0]))
        shell:  "ln -s -r {params.filetolink} {output.r1} && ln -s -r {params.filetolink2} {output.r2}"

else:
    rule simulate_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not rundedup else "DEDUP_FASTQ/{file}_dedup.fastq.gz"
        output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        threads: 1
        params: filetolink = lambda w,input: "{r}".format(r=os.path.abspath(input.r1[0]))
        shell:  "ln -s -r {params.filetolink} {output.r1}"
