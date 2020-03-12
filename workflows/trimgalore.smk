TRIMBIN, TRIMENV = env_bin_from_config2(SAMPLES,config,'TRIMMING')

if paired == 'paired':
    rule trimgalore_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = "TRIMMED_FASTQ/{file}_R1_val_1.fq.gz",
                o2 = "TRIMMED_FASTQ/{file}_R2_val_2.fq.gz"
        log:   "LOGS/{file}_trim.log"
        conda: "snakes/envs/"+TRIMENV+".yaml"
        threads: min(int(MAXTHREAD/2),4) if min(int(MAXTHREAD/2),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
        params: odir=lambda wildcards,output:os.path.dirname(output.o1),
                tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")['OPTIONS'][0].items()),
                trim=TRIMBIN
        shell:  "{params.trim} --cores {threads} --paired --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1} {input.r2} &> {log}"

    rule trimgalore_rename:
        input:  o1 = rules.trimgalore_trim.output.o1,
                o2 = rules.trimgalore_trim.output.o2
        output: r1 = "TRIMMED_FASTQ/{file}_R1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_R2_trimmed.fastq.gz"
        conda: "snakes/envs/"+TRIMENV+".yaml"
        threads: 1
        shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"

else:
    rule trimgalore_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
        log:    "LOGS/{file}_trim.log"
        conda: "snakes/envs/"+TRIMENV+".yaml"
        threads: min(int(MAXTHREAD/2),4) if min(int(MAXTHREAD/2),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
        params: odir=lambda wildcards,output: os.path.dirname(output.o1),
                tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")['OPTIONS'][0].items()),
                trim=TRIMBIN
        shell:  "{params.trim} --cores {threads} --no_report_file --gzip {params.tpara} -o {params.odir} {input.r1} &> {log}"

    rule trimgalore_rename:
        input:  o1 = rules.trimgalore_trim.output.o1
        output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        conda: "snakes/envs/"+TRIMENV+".yaml"
        threads: 1
        shell:  "mv {input.o1} {output.r1}"
