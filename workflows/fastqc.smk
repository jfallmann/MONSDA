QCBIN, QCENV = env_bin_from_config(config, 'PREQC')

if QC_PAIRED == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = lambda wildcards: raw_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = lambda wildcards: dedup_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/{read}_fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = lambda wildcards: trimmed_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/{read}_fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = lambda wildcards: raw_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/fastqc/fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = lambda wildcards: dedup_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = lambda wildcards: trimmed_fastq_input(wildcards)
        output: o1 = report("QC/{combo}/{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

rule qc_mapped:
    input:   r1 = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output:  o1 = report("QC/{combo}/{file}_mapped_sorted_fastqc.zip", category="QC")
    log:     "LOGS/{combo}/{file}/QC/fastqc/fastqc_mapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

rule qc_uniquemapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_uniquemapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

rule qc_dedupmapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_dedupmapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

rule qc_uniquededup:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/fastqc/fastqc_uniquededup.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
