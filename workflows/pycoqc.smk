rule pycoqc_raw:
    input: "FASTQ/{rawfile}.fastq.gz"
    output: report("QC/{rawfile}_pycoqc.zip", category="QC")
    log:    "LOGS/{rawfile}/pycoqc_raw.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input}; do OUT=$(dirname {output});pycoqc --quiet -o $OUT -t {threads} --noextract -f fastq {input} 2> {log};done"

rule pycoqc_trimmed:
    input:  "TRIMMED_FASTQ/{rawfile}_trimmed.fastq.gz",
            "QC/{rawfile}_pycoqc.zip"
    output: report("QC/{rawfile}_trimmed_pycoqc.zip", category="QC")
    log:    "LOGS/{rawfile}/pycoqc_trimmed.log"
    conda:  "../envs/qc.yaml"
    threads: 20
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    shell: "for i in {input[0]}; do OUT=$(dirname {output});pycoqc --quiet -o $OUT -t {threads} --noextract -f fastq {input[0]} 2> {log};done"

rule pyqc_mapped:
    input:  "SORTED_MAPPED/{file}_mapped_sorted.sam.gz"
    output:  report("QC/{file}_mapped_sorted_pycoqc.zip", category="QC")
    log: "LOGS/{file}/pycoqc_mapped.log"
    params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
    conda: "../envs/qc.yaml"
    threads: 20
    shell: "for i in {input}; do OUT=$(dirname {output});pycoqc --quiet -o $OUT -t {threads} --noextract -f sam_mapped {input} 2> {log};done"

rule pyqc_uniquemapped:
    input:  "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam",
            "UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: report("QC/{file}_mapped_sorted_unique_pycoqc.zip", category="QC")
    log: "LOGS/{file}/pycoqc_uniquemapped.log"
    conda: "../envs/qc.yaml"
    threads: 20
    params:  dir=lambda w: expand("QC/{source}",source=source_from_sample(w.file))
#    params: dir=expand("QC/{source}",source=SOURCE)
    shell: "for i in {input[0]}; do OUT=$(dirname {output});pycoqc --quiet -o $OUT -t {threads} --noextract -f bam {input[0]} 2> {log};done"
