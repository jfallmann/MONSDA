#include: "header.smk"

rule qc_raw_paired:
    input: r1 = "FASTQ/{rawfile}_r1.fastq.gz",
           r2 = "FASTQ/{rawfile}_r1.fastq.gz"
    output: report("QC/{rawfile}_qc.zip", category="QC")
    wildcard_constraints:
        file="trimmed{0}"
    log:    "LOGS/{rawfile}/fastqc_raw.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2> {log} && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_trimmed_paired:
    input:  r1 = "TRIMMED_FASTQ/{rawfile}_r1_trimmed.fastq.gz",
            r2 = "TRIMMED_FASTQ/{rawfile}_r2_trimmed.fastq.gz",
            rules.qc_raw_paired.output
    output: report("QC/{rawfile}_trimmed_qc.zip", category="QC")
    log:   "LOGS/{rawfile}/fastqc_trimmed.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile))
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract -f fastq {input.r2} 2> {log} && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_mapped_paired:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  report("QC/{file}_mapped_sorted_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log} && cd $OUT && rename fastqc qc *_fastqc*"

rule qc_uniquemapped_paired:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_qc.zip", category="QC")
    log: "LOGS/{file}/fastqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "OUT=$(dirname {output});fastqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log} && cd $OUT && rename fastqc qc *_fastqc*"

include: "multiqc.smk"
