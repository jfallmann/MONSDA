if paired is 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input: r1 = "FASTQ/{rawfile}_r1.fastq.gz",
               r2 = "FASTQ/{rawfile}_r2.fastq.gz"
        #input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
        #        r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = report("QC/{rawfile}_r1_fastqc.zip",category="QC"),
                o2 = report("QC/{rawfile}_r2_fastqc.zip",category="QC")
        wildcard_constraints:
            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2>> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_trimmed:
        input:  expand(rules.qc_raw_paired.output, rawfile=SAMPLES),
        #input:  raw1 = lambda wildcards: expand(rules.qc_raw.output.o1, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
        #        raw2 = lambda wildcards: expand(rules.qc_raw.output.o2, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                r1 = "TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",
                r2 = "TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz"
        output: o1 = report("QC/{file}_r1_trimmed_fastqc.zip", category="QC"),
                o2 = report("QC/{file}_r2_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2>> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_mapped:
        input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
        output:  report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped:
        input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
                "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
        output: report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    #    params: dir=expand("QC/{source}",source=SOURCE)
        shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log}"

else:
    rule qc_raw:
        input: q1 = "FASTQ/{rawfile}.fastq.gz"
        #input:  r1 = lambda wildcards: "FASTQ/{rawfile}.fastq.gz".format(rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = report("QC/{rawfile}_fastqc.zip", category="QC")
        wildcard_constraints:
            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_trimmed:
        input: expand(rules.qc_raw.output.o1, rawfile=SAMPLES),
        #input:  raw = lambda wildcards: expand(rules.qc_raw.output.o1, rawfile=[x for x in SAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
                q1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz",file=samplecond(SAMPLES,config)),
        output: o1 = report("QC/{file}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_mapped:
        input:   q1 = expand("SORTED_MAPPED/{file}_mapped_sorted.sam.gz",file=samplecond(SAMPLES,config))
        output:  o1 = report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped:
        input:  q1 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",file=samplecond(SAMPLES,config)),
                q2 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai",file=samplecond(SAMPLES,config))
        output: o1 = report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input.q1} 2> {log}"
