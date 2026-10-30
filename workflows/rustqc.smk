QCBIN, QCENV = env_bin_from_config(config, 'POSTQC')

# Map MONSDA strandedness to RustQC strandedness
def rustqc_stranded(stranded):
    if stranded == 'fr':
        return 'forward'
    elif stranded == 'rf':
        return 'reverse'
    else:
        return 'unstranded'

RUSTQC_STRANDED = rustqc_stranded(QC_STRANDED)
RUSTQC_PAIRED = '-p' if QC_PAIRED == 'paired' else ''

rule rustqc_mapped:
    input:   r1 = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output:  o1 = directory("QC/{combo}/{file}_mapped_sorted"),
             js = "QC/{combo}/{file}_mapped_sorted/rustqc_summary.json",
             tmpanno = temp("TMP/QC/{combo}/{file}_ms_anno.gtf")
    log:     "LOGS/{combo}/{file}/QC/rustqc/rustqc_mapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED,
             bins = BINS
    shell: "(gzip -cdfq {params.anno:q} > {output.tmpanno:q} && rustqc rna {input.r1:q} --gtf {output.tmpanno:q} -t {threads} {params.paired} -s {params.stranded:q} --skip-dup-check -j {output.js:q} -o {output.o1:q} {params.qpara} && python {params.bins:q}/Analysis/check_strandedness.py --qc-dir {output.o1:q} --expected {params.stranded:q} --sample {wildcards.file:q} --report {output.o1:q}/strandedness_check.json) 2> {log:q}"

rule rustqc_uniquemapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = directory("QC/{combo}/{file}_mapped_sorted_unique"),
        js = "QC/{combo}/{file}_mapped_sorted_unique/rustqc_summary.json",
        tmpanno = temp("TMP/QC/{combo}/{file}_msu_anno.gtf")
    log:    "LOGS/{combo}/{file}/QC/rustqc/rustqc_uniquemapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED,
             bins = BINS
    shell: "(gzip -cdfq {params.anno:q} > {output.tmpanno:q} && rustqc rna {input.r1:q} --gtf {output.tmpanno:q} -t {threads} {params.paired} -s {params.stranded:q} --skip-dup-check -j {output.js:q} -o {output.o1:q} {params.qpara} && python {params.bins:q}/Analysis/check_strandedness.py --qc-dir {output.o1:q} --expected {params.stranded:q} --sample {wildcards.file:q} --report {output.o1:q}/strandedness_check.json) 2> {log:q}"

rule rustqc_dedupmapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_{type}_dedup.bam"
    output: o1 = directory("QC/{combo}/{file}_mapped_{type}_dedup"),
            js = "QC/{combo}/{file}_mapped_{type}_dedup/rustqc_summary.json",
            tmpanno = temp("TMP/QC/{combo}/{file}_{type}_dedup_anno.gtf")
    log:    "LOGS/{combo}/{file}/QC/rustqc/rustqc_dedupmapped_{type}.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED,
             bins = BINS
    shell: "(gzip -cdfq {params.anno:q} > {output.tmpanno:q} && rustqc rna {input.r1:q} --gtf {output.tmpanno:q} -t {threads} {params.paired} -s {params.stranded:q} --skip-dup-check -j {output.js:q} -o {output.o1:q} {params.qpara} && python {params.bins:q}/Analysis/check_strandedness.py --qc-dir {output.o1:q} --expected {params.stranded:q} --sample {wildcards.file:q} --report {output.o1:q}/strandedness_check.json) 2> {log:q}"
