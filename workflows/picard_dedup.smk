DEDUPBIN, DEDUPENV = env_bin_from_config(config, 'DEDUP')

rule dedupbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam"
    output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
            bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
            met = report("MAPPED/{combo}/{file}_mapped_{type}_dedup_metrics.txt", category="DEDUP"),
            td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
    log:    "LOGS/{combo}/{file}_{type}/dedupbam.log"
    conda:  ""+DEDUPENV+".yaml"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('JAVA', ""),
            dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && {params.dedup} {params.jpara} MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR {output.td} --INPUT {input.bam} --OUTPUT {output.bam} --METRICS_FILE {output.met} {params.dpara} 2> {log} && samtools index {output.bam} 2>> {log}"