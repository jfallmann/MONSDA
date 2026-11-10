DEDUPBIN, DEDUPENV = env_bin_from_config(config, 'DEDUP')

rule dedupbam:
    input:  bam = "MAPPED/{combo}/{file}_mapped_{type}.bam"
    output: bam = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam", category="DEDUP"),
            bai = report("MAPPED/{combo}/{file}_mapped_{type}_dedup.bam.bai", category="DEDUP"),
            met = report("MAPPED/{combo}/{file}_mapped_{type}_dedup_metrics.txt", category="DEDUP"),
            td = temp(directory("TMP/UMIDD/{combo}/{file}_{type}"))
    log:    "LOGS/{combo}/{file}/DEDUP/picard/dedupbam_{type}.log"
    conda:  ""+DEDUPENV+".yaml"
    container: "oras://jfallmann/monsda:"+DEDUPENV+""
    threads: 1
    priority: 0               # This should be done after all mapping is done
    params: jpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('JAVA', ""),
            dpara = lambda wildcards: tool_params(wildcards.file, None, config, "DEDUP", DEDUPENV)['OPTIONS'].get('DEDUP', ""),
            dedup = DEDUPBIN
    shell: "mkdir -p {output.td}; RG={output.td}/{wildcards.file}_rg.bam; : > {log}; if samtools view -H {input.bam} | grep '^@RG' >/dev/null; then ln -sf $(realpath {input.bam}) $RG; else {params.dedup} {params.jpara} AddOrReplaceReadGroups --INPUT {input.bam} --OUTPUT $RG --RGID {wildcards.file} --RGLB MONSDA --RGPL UNKNOWN --RGPU {wildcards.file} --RGSM {wildcards.file} >> {log} 2>&1; fi && {params.dedup} {params.jpara} MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR {output.td} --INPUT $RG --OUTPUT {output.bam} --METRICS_FILE {output.met} {params.dpara} 2>> {log} && rm -f $RG && samtools index {output.bam} 2>> {log} && cp -f {output.met} $(dirname {log})/"