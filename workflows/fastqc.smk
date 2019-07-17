#include: "header.smk"

rule qc_raw:
    input: "FASTQ/{rawfile}.fastq.gz",
    output: report("QC/{rawfile}_fastqc.zip", category="QC")
    wildcard_constraints:
        file="!trimmed"
    log:    "LOGS/{rawfile}/fastqc_raw.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
    shell: "do OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_trimmed:
    input:  "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",
            rules.qc_raw.output
    output: report("QC/{rawfile}_trimmed_fastqc.zip", category="QC")
    log:   "LOGS/{rawfile}/fastqc_trimmed.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input[0]} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_mapped:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  report("QC/{file}_mapped_sorted_fastqc.zip", category="QC")
    log: "LOGS/{file}/fastqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_uniquemapped:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log: "LOGS/{file}/fastqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

include: "multiqc.smk"
