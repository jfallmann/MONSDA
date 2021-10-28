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
        if (filename.indexOf("_dedup.bam") > 0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam"
        else if (filename.indexOf("_dedup.bam.bai") > 0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}_dedup.bam.bai"
        else if (filename.indexOf("log") > 0)    "LOGS/$COMBO$CONDITION/dedupbam.log"
        else null
    }

    input:
    val dummy
    path samples
    
        
    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: log

    script:
    if (PAIRED == 'paired'){
        out=samples.getSimpleName()+"_dedup.bam"
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=ded.log --paired --stdin=$samples --stdout=$out && samtools index $out &>> ded.log
        """
    }
    else{
        out=samples.getSimpleName()+"_dedup.fastq.gz"
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=ded.log --stdin=$samples --stdout=$out && samtools index $out &>> ded.log
        """
    }
    """
    
    """
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

    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true)
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true)
    msamples_ch.combine(usamples_ch)

    collect_dedup(collection.collect())
    dedup(collect_dedup.out.done, msamples_ch)

    emit:
    dedup = dedup.out.bam
}


