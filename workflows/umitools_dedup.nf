DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

DEDUPPARAMS = get_always('umitools_params_DEDUP') ?: ''

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
        if (filename.indexOf("_dedup.bam") > 0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf("_dedup.bam.bai") > 0) "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)           "LOGS/$COMBO$CONDITION/DEDUP/dedupbam.log"
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
    if (PAIRED == 'paired'){        
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=ded.log --paired --stdin=$samples --stdout=$out && samtools index $out &>> ded.log
        """
    }
    else{
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=ded.log --stdin=$samples --stdout=$out && samtools index $out &>> ded.log
        """
    }
}

workflow DEDUPBAM{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/$COMBO"+element+"_R1.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/$COMBO"+element+"_R2.fastq.gz"
        }
        T2SAMPLES.sort()
        dedup_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/$COMBO"+element+".fastq.gz"
        }
        T1SAMPLES.sort()
        dedup_samples_ch = Channel.fromPath(T1SAMPLES)
    }

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
    mindex_ch = Channel.fromPath(MINDEX, followLinks: true)
    uindex_ch = Channel.fromPath(UINDEX, followLinks: true)
    msamples_ch.join(usamples_ch)
    mindex_ch.join(uindex_ch)

    collect_dedup(collection.collect())
    dedup(collect_dedup.out.done, msamples_ch, mindex_ch)

    emit:
    dedup = dedup.out.bam
}


