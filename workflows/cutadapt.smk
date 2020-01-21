TRIMBIN, TRIMENV = env_bin_from_config2(SAMPLES,config,'TRIMMING')

if paired == 'paired':
	rule cutadapt_trim:
	    input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
	            r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
	    output: o1 = "TRIMMED_FASTQ/{file}_r1_val_1.fq.gz",
	            o2 = "TRIMMED_FASTQ/{file}_r2_val_2.fq.gz"
	    log:    "LOGS/{file}_trim.log"
	    conda: "snakes/envs/"+TRIMENV+".yaml"
	    threads: lambda x: min(int(MAXTHREAD/8),4) if min(int(MAXTHREAD/8),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
	    params: odir=lambda wildcards,output:os.path.dirname(output.o1),
	            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
	            trim=TRIMBIN
	    shell:  "{params.trim} -a file:{params.ada} {params.tpara} --cores {threads} -o {output.r1} -p {output.r2} {input.r1} {input.r2} &> {log}"

	rule cutadapt_rename:
	    input:  o1 = rules.cutadapt_trim.output.o1,
	            o2 = rules.cutadapt_trim.output.o2
	    output: r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
	            r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
	    conda: "snakes/envs/"+TRIMENV+".yaml"
	    threads: 1
	    shell:  "mv {input.o1} {output.r1} && mv {input.o2} {output.r2}"

else:
    rule cutadapt_trim:
	    input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
	    output: o1 = "TRIMMED_FASTQ/{file}_trimmed.fq.gz"
	    log:    "LOGS/{file}_trim.log"
	    conda: "snakes/envs/"+TRIMENV+".yaml"
	    threads: lambda x: min(int(MAXTHREAD/8),4) if min(int(MAXTHREAD/2),4) >= 1 else (4 if int(MAXTHREAD) >= 4 else 1)
	    params: odir=lambda wildcards,output: os.path.dirname(output.o1),
	            tpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, "TRIMMING")[0].items()),
	            trim=TRIMBIN,
	    shell:  "{params.trim} -a file:{params.ada} {params.tpara} --cores {threads} -o {output.o1} {input.r1} > {log}"

	rule cutadapt_rename:
	    input:  o1 = rules.cutadapt_trim.output.o1
	    output: r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
	    conda: "snakes/envs/"+TRIMENV+".yaml"
	    threads: 1
	    shell:  "mv {input.o1} {output.r1}"
