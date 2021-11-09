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
        if (filename.endsWith("_dedup.bam"))              "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam"
        else if (filename.indexOf("_dedup.bam.bai") > 0)  "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam.bai"
        else if (filename.indexOf("log") > 0)             "LOGS/$COMBO$CONDITION/DEDUP/dedupbam.log"
        else if (filename.indexOf("metrix.txt") > 0)      "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedupmetrics.txt"
        else null
    }

    input:
    path dummy
    path samples
    path indices
        
    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: log

    script:
    out=samples.getSimpleName()+"_dedup.bam"
    """
    mkdir -p tmp && $DEDUPBIN $JAVAPARAMS MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR tmp --INPUT $samples --OUTPUT $out --METRICS_FILE dedup_metrics.txt $DEDUPPARAMS  &> dedup.log && samtools index $out &>> dedup.log
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
    MINDEX = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted.bam.bai"
    }
    MINDEX.sort()
    
    UINDEX = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted_unique.bam.bai"
    }
    UINDEX.sort()
    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true)
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true)
    msamples_ch.join(usamples_ch)
    mindex_ch = Channel.fromPath(MINDEX, followLinks: true)
    uindex_ch = Channel.fromPath(UINDEX, followLinks: true)
    mindex_ch.join(uindex_ch)
    
    collect_dedup(collection.collect())
    dedup(collect_dedup.out.done, msamples_ch, mindex_ch)

    emit:
    dedup = dedup.out.bam
}


