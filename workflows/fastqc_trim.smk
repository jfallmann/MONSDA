if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/RAW/{rawfile}_{read}_fastqc.zip", category="QC")
        log:    "LOGS/{rawfile}/fastqc_{read}_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}", source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_{read}_trimmed.fastq.gz", file=samplecond(SAMPLES,config), read=["R1", "R2"])
        output: o1 = report("QC/TRIMMED/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_{read}_trimmed.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}", source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=SAMPLES, read=['R1','R2']),
                expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config),read=['R1','R2'])
        output: html = report("QC/Multi/TRIMMED_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/TRIMMED_RAW/{condition}/tmp"),
                lst = "QC/Multi/TRIMMED_RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_trimmed_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/RAW/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{rawfile}/fastqc_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}", source=source_from_sample(w.rawfile,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(SAMPLES[0], None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = expand("TRIMMED_FASTQ/{file}_trimmed.fastq.gz", file=samplecond(SAMPLES,config))
        output: o1 = report("QC/TRIMMED/{file}_trimmed_fastqc.zip", category="QC")
        log:   "LOGS/{file}/fastqc_trimmed.log"
        conda:  "snakes/envs/qc.yaml"
        threads: MAXTHREAD
        params: dir=lambda w: expand("QC/{source}", source=source_from_sample(w.file,config)),
                qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key,val) for (key,val) in tool_params(wildcards.file, None ,config, 'QC')['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=SAMPLES),
               expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES,config))
        output: html = report("QC/Multi/TRIMMED_RAW/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/TRIMMED_RAW/{condition}/tmp"),
                lst = "QC/Multi/TRIMMED_RAW/{condition}/qclist.txt"
        log:    "LOGS/{condition}/multiqc_trimmed_raw.log"
        conda:  "snakes/envs/qc.yaml"
        threads: 1
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o $OUT -l {output.lst} 2> {log}"
