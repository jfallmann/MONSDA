QCBIN, QCENV = env_bin_from_config(config, 'PREQC')
#outdir = 'QC/'+str(QCENV)+'/'
#moutdir = 'QC/Multi/'+str(QCENV)+'/'

#wildcard_constraints:
#    rawfile = '|'.join(list(SAMPLES)),
#    read = "R1|R2"
#    outdir = outdir,
#    moutdir = moutdir

#rule themall:
#    input: expand("QC/Multi/{combo}/{condition}/multiqc_report.html", moutdir=moutdir, condition=str.join(os.sep, conditiononly(SAMPLES[0], config)))

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/rastqc/fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_{read}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/{read}_fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/{read}_fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), read=['R1','R2'], combo=combo),
               expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo),
               expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), read=['R1','R2'], combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_trim_dedup_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc_trim_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=QC/{wildcards.combo}/{wildcards.condition}; mkdir -p $OUT; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -o $OUT $SCAN 2> {log}"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/rastqc/fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule multiqc:
        input: expand(rules.qc_raw.output.o1, rawfile=list(SAMPLES), combo=combo),
               expand(rules.qc_dedup.output.o1, file=samplecond(SAMPLES, config), combo=combo),
               expand(rules.qc_trimmed.output.o1, file=samplecond(SAMPLES, config), combo=combo)
        output: html = report("QC/Multi/{combo}/{condition}/multiqc_trim_dedup_report.html", category="QC")
        log:    "LOGS/{combo}/MULTIQC/multiqc/{condition}_multiqc_trim_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: 1
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('MULTI', "")
        shell:  "OUT=$(dirname {output.html}); SCAN=QC/{wildcards.combo}/{wildcards.condition}; mkdir -p $OUT; export LC_ALL=C.UTF-8; multiqc -f {params.qpara} --exclude picard --exclude gatk -k json -z -o $OUT $SCAN 2> {log}"
