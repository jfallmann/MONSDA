QCBIN, QCENV = env_bin_from_config2(SAMPLES,config,'QC')

wildcard_constraints:
    rawfile = '|'.join(list(SAMPLES)),
    read = "R1|R2"

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/FASTQC/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/QC/{rawfile}_fastqc_{read}_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{file}_{read}_dedup.fastq.gz"
        output: o1 = report("QC/FASTQC/{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/QC/{file}_{read}_fastqc_dedup.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2']),
                expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES,config), read=['R1','R2'])
        output: html = report("QC/Multi/DEDUP_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/DEDUP_RAW/{condition}/tmp"),
                lst = "QC/Multi/DEDUP_RAW/{condition}/qclist.txt"
        log:    "LOGS/QC/{condition}_multiqc_trimmed_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/FASTQC/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/QC/{rawfile}_fastqc_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{file}_dedup.fastq.gz"
        output: o1 = report("QC/FASTQC/{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/QC/{file}_fastqc_dedup.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES)),
		       expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES,config))
        output: html = report("QC/Multi/DEDUP_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/DEDUP_RAW/{condition}/tmp"),
                lst = "QC/Multi/DEDUP_RAW/{condition}/qclist.txt"
        log:    "LOGS/QC/{condition}_multiqc_trimmed_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
