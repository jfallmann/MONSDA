if rundedup:
    if paired == 'paired':
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2']),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config),read=['R1','R2']),
                    expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES,config),read=['R1','R2']),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES,config)),
                    expand(rules.dedupuniqbam.output.bam, file=samplecond(SAMPLES,config))
            output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{condition}/tmp"),
                    lst = "QC/Multi/{condition}/qclist.txt"
            log:    "LOGS/{condition}/multiqc.log"
            conda:  "nextsnakes/envs/qc.yaml"
            threads: 1
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

    else:
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES)),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES,config)),
                    expand(rules.dedupbam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.dedupuniqbam.output.bam, file=samplecond(SAMPLES,config))
            output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{condition}/tmp"),
                    lst = "QC/Multi/{condition}/qclist.txt"
            log:    "LOGS/{condition}/multiqc.log"
            conda:  "nextsnakes/envs/qc.yaml"
            threads: 1
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    if paired == 'paired':
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2']),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config),read=['R1','R2']),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES,config))
            output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{condition}/tmp"),
                    lst = "QC/Multi/{condition}/qclist.txt"
            log:    "LOGS/{condition}/multiqc.log"
            conda:  "nextsnakes/envs/qc.yaml"
            threads: 1
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

    else:
        rule multiqc:
            input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES)),
                    expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_mapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.qc_uniquemapped.output.o1, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bam.output.bam, file=samplecond(SAMPLES,config)),
                    expand(rules.sam2bamuniq.output.uniqbam, file=samplecond(SAMPLES,config))
            output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                    tmp = temp("QC/Multi/{condition}/tmp"),
                    lst = "QC/Multi/{condition}/qclist.txt"
            log:    "LOGS/{condition}/multiqc.log"
            conda:  "nextsnakes/envs/qc.yaml"
            threads: 1
            shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
