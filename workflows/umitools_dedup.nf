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

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf("_dedup.bam.bai") > 0) "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)           "LOGS/$COMBO$CONDITION/DEDUP/dedupbam.log"
        else null
    }

    input:
    path bams
    path bais
      
    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: log

    script:
    out=bams.getSimpleName()+"_dedup.bam"
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
 
    bams = collection.filter(~/.bam/)
    bais = collection.filter(~/.bai/)
    dedup(bams, bais)

    emit:
    dedup = dedup.out.bam
}


