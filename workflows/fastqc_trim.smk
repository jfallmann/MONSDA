QCBIN, QCENV = env_bin_from_config3( config, 'QC')
#outdir = 'QC/'+str(QCENV)+'/'
#moutdir = 'QC/Multi/'+str(QCENV)+'/'

#wildcard_constraints:
#    rawfile = '|'.join(list(SAMPLES)),
#    read = "R1|R2"
#    outdir = outdir,
#    moutdir = moutdir

#rule themall:
#    input: expand("{moutdir}TRIMMED_RAW/{condition}/multiqc_report.html", moutdir=moutdir, condition=str.join(os.se#p, conditiononly(SAMPLES[0], config)))

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{combo}/{rawfile}_fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}_{read}_fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input:  expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
                expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                lst = "QC/Multi/{combo}/{condition}/qclist_trimmed_raw.txt"
        log:    "LOGS/{combo}/{condition}_multiqc_trimmed_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: 1
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[1], None, config, 'QC', QCENV)['OPTIONS'][1].items())
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/QC/{combo}/{rawfile}_fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/QC/{combo}/{file}_fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'][0].items())
        shell: "OUT=$(dirname {output.o1});fastqc --quiet -o $OUT -t {threads} --noextract {params.qpara} -f fastq {input.r1} 2> {log}"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
               expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_report.html", category="QC"),
                tmp = temp("QC/Multi/{combo}/{condition}/tmp"),
                lst = "QC/Multi/{combo}/{condition}/qclist_trimmed_raw.txt"
        log:    "LOGS/{combo}/{condition}_multiqc_trimmed_raw.log"
        conda:  ""+QCENV+".yaml"
        threads: 1
        params:  qpara = lambda wildcards: ' '.join("{!s} {!s}".format(key, val) for (key, val) in tool_params(SAMPLES[1], None, config, 'QC', QCENV)['OPTIONS'][1].items())
        shell:  "OUT=$(dirname {output.html}); for i in {input};do echo $(dirname \"${{i}}\") >> {output.tmp};done; cat {output.tmp} |sort -u > {output.lst};export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -s -o $OUT -l {output.lst} 2> {log}"
