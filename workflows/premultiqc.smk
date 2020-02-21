rule qcall:
    input: expand("QC/Multi/{condition}/multiqc_report.html",condition=os.path.join(*samplecond(SAMPLES,config)[0].split(os.sep)[:-1]))

if paired == 'paired':
    rule multiqc:
        input: expand("QC/{rawfile}_{read}_fastqc.zip", rawfile=SAMPLES, read=['r1','r2']),
        output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{condition}/tmp"),
                lst = "QC/Multi/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    rule multiqc:
        input: expand("QC/{rawfile}_fastqc.zip", rawfile=SAMPLES),
        output: html = report("QC/Multi/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{condition}/tmp"),
                lst = "QC/Multi/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
