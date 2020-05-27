if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}_{read}.fastq.gz", rawfile=list(SAMPLES), read=['R1','R2'])
        output: o1 = report(expand("QC/FASTQC/{rawfile}_{read}_fastqc.zip", rawfile=list(SAMPLES), read=['R1','R2']), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_{read}_raw.log", rawfile=list(SAMPLES), read=['R1','R2'])
        conda:  "nextsnakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2> {log};done"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'])
        output: html = report("QC/Multi/RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/RAW/{condition}/tmp"),
                lst = "QC/Multi/RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_raw.log"
        conda:  "nextsnakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = expand("FASTQ/{rawfile}.fastq.gz", rawfile=list(SAMPLES))
        output: o1 = report(expand("QC/FASTQC/{rawfile}_fastqc.zip", rawfile=list(SAMPLES)), category="QC")
        log:    expand("LOGS/{rawfile}/fastqc_raw.log", rawfile=list(SAMPLES))
        conda:  "nextsnakes/envs/qc.yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "arr=({input.r1}); orr=({output.o1});alen=${{#arr[@]}}; for i in \"${{!arr[@]}}\";do OUT=$(dirname ${{orr[$i]}}); fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq  ${{arr[$i]}} 2>> {log};done"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES))
        output: html = report("QC/Multi/RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/RAW/{condition}/tmp"),
                lst = "QC/Multi/RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_raw.log"
        conda:  "nextsnakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
