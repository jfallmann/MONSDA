#include: "header.smk"

rule qc_raw:
    input: q1 = "FASTQ/{rawfile}.fastq.gz",
    output: o1 = report("QC/{rawfile}_fastqc.zip", category="QC")
    wildcard_constraints:
        file="!trimmed"
    log:    "LOGS/{rawfile}/fastqc_raw.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_trimmed:
    input: expand(rules.qc_raw.output.o1, rawfile=SAMPLES),
           q1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
    output: o1 = report("QC/{file}_trimmed_fastqc.zip", category="QC")
    log:   "LOGS/{file}/fastqc_trimmed.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_mapped:
    input:  q1 = "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  o1 = report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
    log: "LOGS/{file}/fastqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_uniquemapped:
    input:  q1 = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            q2 = "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log: "LOGS/{file}/fastqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input.q1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

include: "multiqc.smk"
