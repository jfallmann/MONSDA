TRIMBIN, TRIMENV = env_bin_from_config2(MAPSAMPLES, config,'TRIMMING')
outdir = 'TRIMMED_FASTQ/'+str(TRIMENV)+'/'

wildcard_constraints:
    rawfile = '|'.join(list(SAMPLES)),
    read = "R1|R2",
    outdir = outdir

#rule trimthemall:
#    input: expand("{outdir}{file}_{read}_trimmed.fastq.gz", outdir=outdir, file=samplecond(SAMPLES, config), read=["R1","R2"]) if paired == \'paired\' else expand("{outdir}{file}_trimmed.fastq.gz", outdir=outdir, file=samplecond(SAMPLES, config))

if paired = 'paired':
    rule bbduk_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_R1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R1_dedup.fastq.gz",
                r2 = lambda wildcards: "FASTQ/{rawfile}_R2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_R2_dedup.fastq.gz"
        output: o1 = "{outdir}{file}_R1_val_1.fq.gz",
                o2 = "{outdir}{file}_R2_val_2.fq.gz"
        log:    "LOGS/{combo}/{file}_trim.log"
        conda: ""+TRIMENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards, output:os.path.dirname(output.o1),
                tpara = lambda wildcards: tool_params(wildcards.file, None, config, "TRIMMING", TRIMENV).get('TRIM', ""),
                trim=TRIMBIN
        shell:  "{params.trim} t={threads} in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} ref={params.ada} {params.tpara}"

    rule bbduk_rename:
        input:  o1 = rules.bbduk_trim.output.o1,
                o2 = rules.bbduk_trim.output.o2
        output: r1 = "{outdir}{file}_R1_trimmed.fastq.gz",
                r2 = "{outdir}{file}_R2_trimmed.fastq.gz"
        conda: ""+TRIMENV+".yaml"
        threads: 1
        shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"
else:
    rule bbduk_trim:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]) if not prededup else "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: o1 = "{outdir}{file}_trimmed.fq.gz"
        log:    "LOGS/{combo}/{file}_trim.log"
        conda: ""+TRIMENV+".yaml"
        threads: MAXTHREAD
        params: odir=lambda wildcards, output: os.path.dirname(output.o1),
                tpara = lambda wildcards: tool_params(wildcards.file, None, config, "TRIMMING", TRIMENV).get('TRIM',""),
                trim=TRIMBIN,
        shell:  "{params.trim} t={threads} in={input.r1} out={output.o1} ref={params.ada} {params.tpara}"

    rule bbduk_rename:
        input:  o1 = rules.bbduk_trim.output.o1
        output: r1 = "{outdir}{file}_trimmed.fastq.gz"
        conda: ""+TRIMENV+".yaml"
        threads: 1
        shell:  "mv {input.o1} {output.r1}"
