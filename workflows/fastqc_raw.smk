rule qcall:
    input: expand("QC/Multi/RAW/{condition}/multiqc_report.html",condition=os.path.join(samplecond(SAMPLES,config)[0]))
    #input: expand("QC/Multi/RAW/{condition}/multiqc_report.html",condition=os.path.join(*samplecond(SAMPLES,config)[0].split(os.sep)[:-1]))

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input: r1 = "FASTQ/{rawfile}_r1.fastq.gz",
               r2 = "FASTQ/{rawfile}_r2.fastq.gz"
        output: o1 = report("QC/{rawfile}_r1_fastqc.zip",category="QC"),
                o2 = report("QC/{rawfile}_r2_fastqc.zip",category="QC")
#        wildcard_constraints:
#            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log} && fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r2} 2>> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule multiqc:
        input: expand("QC/{rawfile}_{read}_fastqc.zip", rawfile=SAMPLES, read=['r1','r2']),
        output: html = report("QC/Multi/RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/RAW/{condition}/tmp"),
                lst = "QC/Multi/RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"


else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{rawfile}_fastqc.zip", category="QC")
#        wildcard_constraints:
#            rawfile="!trimmed"
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}",source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"#" && cd $OUT && rename fastqc qc *_fastqc*"

    rule multiqc:
        input: expand("QC/{rawfile}_fastqc.zip", rawfile=SAMPLES),
        output: html = report("QC/Multi/RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/RAW/{condition}/tmp"),
                lst = "QC/Multi/RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
