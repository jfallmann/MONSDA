QCBIN, QCENV = env_bin_from_config3(SAMPLES, config, 'QC')
outdir = 'QC/'+str(QCENV)+'/'

wildcard_constraints:
    rawfile = '|'.join(list(SAMPLES)),
    read = "R1|R2"
    outdir = outdir
    moutdir = moutdir

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("{outdir}{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{outdir}{rawfile}_fastqc_{read}_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{file}_{read}_dedup.fastq.gz"
        output: o1 = report("{outdir}{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{outdir}{file}_{read}_fastqc_dedup.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("{outdir}{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{outdir}{file}_{read}_fastqc_trimmed.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("{outdir}{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{outdir}{rawfile}_fastqc_raw.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{file}_dedup.fastq.gz"
        output: o1 = report("{outdir}{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{outdir}{file}_fastqc_dedup.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{file}_trimmed.fastq.gz"
        output: o1 = report("{outdir}{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{outdir}{file}_fastqc_trimmed.log"
        conda:  "nextsnakes/envs/"+QCENV+".yaml"
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

rule qc_mapped:
    input:   r1 = "MAPPED/{file}_mapped_sorted.bam"
    output:  o1 = report("QC/FASTQC/{file}_mapped_sorted_fastqc.zip", category="QC")
    log:     "LOGS/QC/{file}_fastqc_mapped.log"
    conda:  "nextsnakes/envs/"+QCENV+".yaml"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f sam_mapped {input.r1} 2> {log}"

rule qc_uniquemapped:
    input:  r1 = "MAPPED/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = report("QC/FASTQC/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log:    "LOGS/QC/{file}_fastqc_uniquemapped.log"
    conda:  "nextsnakes/envs/"+QCENV+".yaml"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

rule qc_dedupmapped:
    input:  r1 = "MAPPED/{file}_mapped_sorted_dedup.bam",
            r2 = "MAPPED/{file}_mapped_sorted_dedup.bam.bai"
    output: o1 = report("QC/FASTQC/{file}_mapped_sorted_dedup_fastqc.zip", category="QC")
    log:    "LOGS/QC/{file}_fastqc_dedupmapped.log"
    conda:  "nextsnakes/envs/"+QCENV+".yaml"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"

rule qc_uniquededup:
    input:  r1 = "MAPPED/{file}_mapped_sorted_unique_dedup.bam",
            r2 = "MAPPED/{file}_mapped_sorted_unique_dedup.bam.bai"
    output: o1 = report("QC/FASTQC/{file}_mapped_sorted_unique_dedup_fastqc.zip", category="QC")
    log:    "LOGS/QC/{file}_fastqc_uniquededup.log"
    conda:  "nextsnakes/envs/"+QCENV+".yaml"
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
    shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f bam {input.r1} 2> {log}"
