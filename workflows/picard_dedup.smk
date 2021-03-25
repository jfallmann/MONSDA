DEDUPBIN, DEDUPENV = env_bin_from_config3(config, 'DEDUP')

rule dedupbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_sorted.bam"
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted_dedup.bam", category="DEDUP"),
            td = temp(directory("TMP/UMIDD/{combo}/{file}"))
    log:    "LOGS/{combo}/{file}/dedupbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items()),
            dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && java {params.jpara} -jar picard.jar MarkDuplicates {params.dedup} --REMOVE_DUPLICATES --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"

rule dedupuniqbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_sorted_unique.bam",
            check = rules.dedupbam.output.bam
    output: bam = report("MAPPED/{combo}/{file}_mapped_sorted_unique_dedup.bam", category="DEDUP"),
            td = temp(directory("TMP/UMIDU/{combo}/{file}"))
    log:    "LOGS/{combo}/{file}/dedupuniqbam.log"
    conda:  "nextsnakes/envs/"+DEDUPENV+".yaml"
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][0].items()),
            dpara = lambda wildcards: ' '.join("{!s}={!s}".format(key, val) for (key, val) in tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'][1].items()),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td} && java {params.jpara} -jar picard.jar MarkDuplicates {params.dedup} --REMOVE_DUPLICATES TRUE --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"
