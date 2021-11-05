DEDUPBIN, DEDUPENV = env_bin_from_config3(config, 'DEDUP')

rule dedupbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted_dedup.bam", category="DEDUP"),
            td = temp(directory("TMP/UMIDD/{combo}/{file}"))
    log:    "LOGS/{combo}/{file}/dedupbam.log"
    conda:  ""+DEDUPENV+".yaml"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('JAVA', ""),
            dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && {params.dedup} {params.jpara} --REMOVE_DUPLICATES --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"

rule dedupuniqbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            check = rules.dedupbam.output.bam
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam", category="DEDUP"),
            td = temp(directory("TMP/UMIDU/{combo}/{file}"))
    log:    "LOGS/{combo}/{file}/dedupuniqbam.log"
    conda:  ""+DEDUPENV+".yaml"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('JAVA', ""),
            dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && {params.dedup} {params.jpara} --REMOVE_DUPLICATES TRUE --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"
