QCBIN, QCENV = env_bin_from_config(config, 'PREQC')

if paired == 'paired':
    log.info('Running paired mode QC')
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}_{read}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_{read}_fastqc.zip")
        log:    "LOGS/{combo}/{rawfile}/QC/rastqc/fastqc_{read}_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_{read}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/{read}_fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_{read}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_{read}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/{read}_fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

else:
    rule qc_raw:
        input:  r1 = "FASTQ/{rawfile}.fastq.gz"
        output: o1 = report("QC/{combo}/{rawfile}_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{rawfile}/QC/rastqc/fastqc_raw.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_dedup:
        input:  r1 = "DEDUP_FASTQ/{combo}/{file}_dedup.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_dedup_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_dedup.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

    rule qc_trimmed:
        input:  r1 = "TRIMMED_FASTQ/{combo}/{file}_trimmed.fastq.gz"
        output: o1 = report("QC/{combo}/{file}_trimmed_fastqc.zip", category="QC")
        log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_trimmed.log"
        conda:  ""+QCENV+".yaml"
        container: "oras://jfallmann/monsda:"+QCENV+""
        threads: MAXTHREAD
        priority: 10
        params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
        shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

rule qc_mapped:
    input:   r1 = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output:  o1 = report("QC/{combo}/{file}_mapped_sorted_fastqc.zip", category="QC")
    log:     "LOGS/{combo}/{file}/QC/rastqc/fastqc_mapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

rule qc_uniquemapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_uniquemapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

rule qc_dedupmapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_dedupmapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"

rule qc_uniquededup:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam.bai"
    output: o1 = report("QC/{combo}/{file}_mapped_sorted_unique_dedup_fastqc.zip", category="QC")
    log:    "LOGS/{combo}/{file}/QC/rastqc/fastqc_uniquededup.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', "")
    shell: "OUT=$(dirname {output.o1});rastqc --quiet --extract -o $OUT -t {threads} {params.qpara} {input.r1} 2> {log}; FILE=$(basename {input.r1}); BASE=${{FILE%.gz}}; BASE=${{BASE%.*}}_fastqc; sed -i -E '/^>>.*\t(PASS|WARN|FAIL)$/ s/\tPASS$/\tpass/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tWARN$/\twarn/; /^>>.*\t(PASS|WARN|FAIL)$/ s/\tFAIL$/\tfail/' \"$OUT/$BASE/fastqc_data.txt\"; rm \"$OUT/$BASE.zip\"; (cd \"$OUT\" && python -m zipfile -c \"$BASE.zip\" \"$BASE\")"
