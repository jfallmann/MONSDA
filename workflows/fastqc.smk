#include: "header.smk"

rule qc_raw:
    input: expand("FASTQ/{file}.fastq.gz",file=SAMPLES)
#    input: "FASTQ/{file}.fastq.gz"
    output: report(expand("QC/{file}_qc.zip",file=SAMPLES), category="QC")
    wildcard_constraints:
        file="trimmed{0}"
    log:    expand("LOGS/{file}/fastqc_raw.log", file=SAMPLES)
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_trimmed:
    input:  expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz",file=SAMPLES),
#    input:  "TRIMMED_FASTQ/{file}_trimmed.fastq.gz",
            rules.qc_raw.output
    output: report(expand("QC/{file}_trimmed_qc.zip",file=SAMPLES), category="QC")
#    wildcard_constraints:
#        file="trimmed+"
    log:    expand("LOGS/{file}/fastqc_trimmed.log",file=SAMPLES)
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input[0]} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_mapped:
#    input:  rules.mapping.output
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  report("QC/{file}_mapped_sorted_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "for i in {input}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_uniquemapped:
#    input: rules.uniqsam.output
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "for i in {input[0]}; do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log};done && cd $OUT && rename fastqc qc *_fastqc*"
