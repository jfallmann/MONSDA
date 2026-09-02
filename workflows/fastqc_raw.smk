QCBIN, QCENV = env_bin_from_config(config, 'PREQC')

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=QC/{wildcards.combo}/{wildcards.condition}; mkdir -p $OUT; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT $SCAN 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=QC/{wildcards.combo}/{wildcards.condition}; mkdir -p $OUT; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT $SCAN 2> {log}"
