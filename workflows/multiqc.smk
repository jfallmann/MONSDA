#include: "header.smk"

rule multiqc:
    input: expand("QC/{qcfile}_qc.zip", qcfile=SAMPLES),
           expand("QC/{qcfile}_trimmed_qc.zip", qcfile=SAMPLES),
           expand("QC/{file}_mapped_sorted_qc.zip", file=samplecond(SAMPLES,config)),
           expand("QC/{file}_mapped_sorted_unique_qc.zip", file=samplecond(SAMPLES,config)),
           expand("SORTED_MAPPED/{file}_mapped_sorted.bam", file=samplecond(SAMPLES,config)),
           expand("UNIQUE_MAPPED/{file}_mapped_sorted_unique.bam", file=samplecond(SAMPLES,config))
    output: report("QC/Multi/multiqc_report.html", category="QC"),
            temp("QC/Multi/tmp"),
            "QC/Multi/qclist.txt"
    log:    "LOGS/multiqc.log"
    conda:  "../envs/qc.yaml"
    shell:  "OUT=$(dirname {output[0]}); for i in {input};do echo $(dirname \"${{i}}\") >> {output[1]};done; cat {output[1]} |sort -u > {output[2]};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8;multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT <(cat ${output[2]}) 2> {log}"
