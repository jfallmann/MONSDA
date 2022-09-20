QCBIN, QCENV = env_bin_from_config(config, 'QC')

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC{combo}{rawfile}_{read}_fastqc.zip")
        log:    "LOGS{combo}{rawfile}_fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo)
        output: html = report("QC/Multi{combo}{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi{combo}{condition}/tmp"),
                lst = "QC/Multi{combo}{condition}/qclist_raw.txt"
        log:    "LOGS{combo}{condition}_multiqc_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC{combo}{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS{combo}{rawfile}_fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo)
        output: html = report("QC/Multi{combo}{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi{combo}{condition}/tmp"),
                lst = "QC/Multi{combo}{condition}/qclist_raw.txt"
        log:    "LOGS{combo}{condition}_multiqc_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
