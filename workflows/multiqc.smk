#include: "header.smk"

rule multiqc:
    input:  rules.qc_raw.output,
            rules.qc_trimmed.output,
            rules.qc_mapped.output,
            rules.qc_uniquemapped.output,
            rules.sam2bam.output,
            rules.sam2bamuniq.output
    output: report("QC/Multi/DONE", category="QC"),
            temp("QC/Multi/qclist.txt")
    log:    "LOGS/multiqc.log"
    conda:  "../envs/qc.yaml"
    shell:  "OUT=$(dirname {output}); for i in {input};do echo $i > {output[1]}; export LC_ALL=C.UTF-8;multiqc --file-list {output[1]} -k json -z -o $OUT 2> {log} && touch QC/Multi/DONE"
