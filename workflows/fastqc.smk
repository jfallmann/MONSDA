QCSAMPLES=list(set(samples(config)))

if any('_r1' in s for s in QCSAMPLES) and any('_r2' in s for s in QCSAMPLES):
    rule qc_raw_paired:
        input: r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in QCAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
               r2 = lambda wildcards: "FASTQ/{rawfile}_r2.fastq.gz".format(rawfile=[x for x in QCSAMPLES if x.split(os.sep)[-1] in wildcards.file][0])
        output: o1 = report("QC/{rawfile}_r1_fastqc.zip",category="QC"),
                o2 = report("QC/{rawfile}_r2_fastqc.zip",category="QC")
        wildcard_constraints:
            file="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "../qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2>> {log}"

    rule qc_trimmed_paired:
        input:  expand(rules.qc_raw_paired.output, rawfile=QCSAMPLES),
                r1 = expand("TRIMMED_FASTQ/{file}_r1_trimmed.fastq.gz",file=samplecond(QCSAMPLES,config)),
                r2 = expand("TRIMMED_FASTQ/{file}_r2_trimmed.fastq.gz",file=samplecond(QCSAMPLES,config)))
        output: o1 = report("QC/{file}_r1_trimmed_fastqc.zip", category="QC"),
                o2 = report("QC/{file}_r2_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2>> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_mapped_paired:
        input:  expand("SORTED_MAPPED/{file}_mapped_sorted.sam.gz",file=samplecond(QCSAMPLES,config)),
        output:  report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped_paired:
        input:  expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",file=samplecond(QCSAMPLES,config),
                expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai",file=samplecond(QCSAMPLES,config)
        output: report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    #    params: dir=expand("QC/{source}",source=SOURCE)
        shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"
    include: "multiqc_paired.smk"

else:
    rule qc_raw:
        input:  r1 = lambda wildcards: "FASTQ/{rawfile}_r1.fastq.gz".format(rawfile=[x for x in QCAMPLES if x.split(os.sep)[-1] in wildcards.file][0]),
        output: o1 = report("QC/{rawfile}_fastqc.zip", category="QC")
        wildcard_constraints:
            file="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_trimmed:
        input:  expand(rules.qc_raw.output.o1, rawfile=QCSAMPLES),
                q1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz",file=samplecond(QCSAMPLES,config)),
        output: o1 = report("QC/{file}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "../envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_mapped:
        input:   q1 = expand("SORTED_MAPPED/{file}_mapped_sorted.sam.gz",file=samplecond(QCSAMPLES,config))
        output:  o1 = report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_mapped.log"
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule qc_uniquemapped:
        input:  q1 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",file=samplecond(QCSAMPLES,config),
                q2 = expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"file=samplecond(QCSAMPLES,config)
        output: o1 = report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
        log: "LOGS/{file}/fastqc_uniquemapped.log"
        conda: "../envs/qc.yaml"
        threads: MAXTHREAD
        params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    #    params: dir=expand("QC/{source}",source=SOURCE)
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    include: "multiqc.smk"
