DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')
DEDUPPARAMS = get_always('picard_params_DEDUP') ?: ''
JAVAPARAMS = get_always('picard_params_JAVA') ?: ''

process collect_multi{
    input:
    path check
    val checker

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}


process multiqc{
    conda "$DEDUPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "DEDUP/Multi/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "DEDUP/Multi/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    path others
    path samples
    path tsamples
    path msamples
    path usamples
    path logs

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -s .
    """
}

workflow MULTIDEDUP{
    take:
    otherqcs
    maplogs

    main:

    //SAMPLE CHANNELS
    RSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../DEDUP/$COMBO"+element+"_fastqc.zip"
    }
    RSAMPLES.sort()
    samples_ch = Channel.fromPath(RSAMPLES, followLinks: true)

    TSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../DEDUP/$COMBO"+element+"_trimmed_fastqc.zip"
    }
    TSAMPLES.sort()
    tsamples_ch = Channel.fromPath(TSAMPLES, followLinks: true)

    MSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../DEDUP/$COMBO"+element+"_mapped_sorted_fastqc.zip"
    }
    MSAMPLES.sort()

    USAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../DEDUP/$COMBO"+element+"_mapped_sorted_unique_fastqc.zip"
    }
    USAMPLES.sort()

    MAPLOG = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+".log"
    }
    MAPLOG.sort()

    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true)
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true)
    logs_ch = Channel.fromPath(MAPLOG, followLinks: true)

    collect_multi(otherqcs.collect(), maplogs.collect())
    multiqc(collect_multi.out.done.collect(), samples_ch, tsamples_ch, msamples_ch, usamples_ch, logs_ch)

    emit:
    mqcres = multiqc.out.multiqc_results
}


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
    shell: "mkdir -p {output.td} && java {params.jpara} -jar picard.jar MarkDuplicates {params.dedup} --REMOVE_DUPLICATES --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"

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
    shell: "mkdir -p {output.td} && java {params.jpara} -jar picard.jar MarkDuplicates {params.dedup} --REMOVE_DUPLICATES TRUE --ASSUME_SORTED --TMP_DIR={output.td} INPUT={input.bam} OUTPUT={output.bam} 2>> {log}"
