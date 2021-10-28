DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')
DEDUPPARAMS = get_always('picard_params_DEDUP') ?: ''
JAVAPARAMS = get_always('picard_params_JAVA') ?: ''

process collect_dedup{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}


process dedup{
    conda "$DEDUPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("_dedup.bam") > 0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam",
        else if (filename.indexOf("_dedup.bam.bai") > 0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam.bai"
        else if (filename.indexOf("log") > 0)    "LOGS/$COMBO$CONDITION/dedupbam.log"
        else null
    }

    input:
    val dummy
    path samples
        
    output:
    path "*.bam", emit: bam
    path "*.log", emit: log

    script:
    """
    mkdir -p tmp && java $JAVAPARAMS -jar picard.jar MarkDuplicates $DEDUPPARAMS --REMOVE_DUPLICATES --ASSUME_SORTED --TMP_DIR=tmp INPUT=$samples OUTPUT=$samples"_dedup.bam" &> dedup.log && samtools index $out &>> dedup.log
    """
}

workflow DEDUPBAM{
    take:
    collection

    main:

    //SAMPLE CHANNELS
    MSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted.bam"
    }
    MSAMPLES.sort()

    USAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted_unique.bam"
    }
    USAMPLES.sort()

    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true)
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true)
    msamples_ch.combine(usamples_ch)

    collect_dedup(collection.collect())
    dedup(collect_dedup.out.done.collect(), msamples_ch)

    emit:
    dedup = dedup.out.bam
}


