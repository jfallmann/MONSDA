if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input: r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{rawfile}_{read}_fastqc.zip",category="QC")
#        wildcard_constraints:
#            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_{read}_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  expand(rules.qc_raw.output, rawfile=SAMPLES, read=['R1','R2']),
                r1 = "TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("QC/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_{read}_trimmed.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_mapped:
        input:  r1 = "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
        output: o1 = report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        conda: "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f sam_mapped {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped:
        input:  r1 = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
                r2 = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
        output: o1 = report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                 qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
    #    params: dir=expand("QC/{source}",source=SOURCE)
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{rawfile}_fastqc.zip", category="QC")
#        wildcard_constraints:
#            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_trimmed:
        input:  expand(rules.qc_raw.output.o1, rawfile=SAMPLES),
                r1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz",file=samplecond(SAMPLES,config)),
        output: o1 = report("QC/{file}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_mapped:
        input:   r1 = expand("SORTED_MAPPED/{file}_mapped_sorted.sam.gz",file=samplecond(SAMPLES,config))
        output:  o1 = report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        conda: "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f sam_mapped {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped:
        input:  r1 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",file=samplecond(SAMPLES,config)),
                r2 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai",file=samplecond(SAMPLES,config))
        output: o1 = report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file,config)),
                 qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
