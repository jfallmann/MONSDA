QCBIN, QCENV = env_bin_from_config(config, 'POSTQC')

# Map MONSDA strandedness to RustQC strandedness
def rustqc_stranded(stranded):
    if stranded == 'fr':
        return 'forward'
    elif stranded == 'rf':
        return 'reverse'
    else:
        return 'unstranded'

RUSTQC_STRANDED = rustqc_stranded(stranded)
RUSTQC_PAIRED = '-p' if paired == 'paired' else ''

rule rustqc_mapped:
    input:   r1 = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output:  o1 = directory("QC/{combo}/{file}_mapped_sorted"),
             js = "QC/{combo}/{file}_mapped_sorted/rustqc_summary.json"
    log:     "LOGS/{combo}/{file}_rustqc_mapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED
    shell: "rustqc rna {input.r1} --gtf {params.anno} -t {threads} {params.paired} -s {params.stranded} --skip-dup-check -j {output.js} -o {output.o1} {params.qpara} 2> {log}"

rule rustqc_uniquemapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam.bai"
    output: o1 = directory("QC/{combo}/{file}_mapped_sorted_unique"),
        js = "QC/{combo}/{file}_mapped_sorted_unique/rustqc_summary.json"
    log:    "LOGS/{combo}/{file}_rustqc_uniquemapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED
    shell: "rustqc rna {input.r1} --gtf {params.anno} -t {threads} {params.paired} -s {params.stranded} --skip-dup-check -j {output.js} -o {output.o1} {params.qpara} 2> {log}"

rule rustqc_dedupmapped:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_dedup.bam.bai"
    output: o1 = directory("QC/{combo}/{file}_mapped_sorted_dedup"),
        js = "QC/{combo}/{file}_mapped_sorted_dedup/rustqc_summary.json"
    log:    "LOGS/{combo}/{file}_rustqc_dedupmapped.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED
    shell: "rustqc rna {input.r1} --gtf {params.anno} -t {threads} {params.paired} -s {params.stranded} -j {output.js} -o {output.o1} {params.qpara} 2> {log}"

rule rustqc_uniquededup:
    input:  r1 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam",
            r2 = "MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam.bai"
    output: o1 = directory("QC/{combo}/{file}_mapped_sorted_unique_dedup"),
        js = "QC/{combo}/{file}_mapped_sorted_unique_dedup/rustqc_summary.json"
    log:    "LOGS/{combo}/{file}_rustqc_uniquededup.log"
    conda:  ""+QCENV+".yaml"
    container: "oras://jfallmann/monsda:"+QCENV+""
    threads: MAXTHREAD
    params:  qpara = lambda wildcards: tool_params(SAMPLES[0], None, config, 'QC', QCENV)['OPTIONS'].get('QC', ""),
             anno = ANNOTATION,
             paired = RUSTQC_PAIRED,
             stranded = RUSTQC_STRANDED
    shell: "rustqc rna {input.r1} --gtf {params.anno} -t {threads} {params.paired} -s {params.stranded} -j {output.js} -o {output.o1} {params.qpara} 2> {log}"
