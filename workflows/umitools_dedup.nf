DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

DEDUPPARAMS = get_always('umitools_params_DEDUP') ?: ''

process dedup_bam{
    conda "$DEDUPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf("_dedup.bam.bai") > 0) "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf("dedup.log") > 0)           "LOGS/$COMBO$CONDITION/DEDUP/dedupbam.log"
        else null
    }

    input:
    path todedup
    path bami
        
    output:
    path "*_dedup.bam", emit: bam
    path "*_dedup.bam.bai", emit: bai
    path "dedup.log", emit: logs

    script:
    bams = todedup[0]
    bais = todedup[1]
    outf = bams.getSimpleName()+"_dedup.bam"
    if (PAIRED == 'paired'){        
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=dedup.log --paired --stdin=$bams --stdout=$outf && samtools index $outf &> tmp.log && cat tmp.log >> dedup.log
        """
    }
    else{
        """
            mkdir tmp && $DEDUPBIN dedup $DEDUPPARAMS --temp-dir tmp --log=dedup.log --stdin=$bams --stdout=$outf && samtools index $outf &> tmp.log && cat tmp.log >> dedup.log
        """
    }
}

workflow DEDUPBAM{
    take: 
    map
    mapi
    mapu
    mapui

    main:
    dedup_bam(map.concat(mapu), mapi.concat(mapui))

    emit:
    dedup = dedup_bam.out.bam
    dedupbai = dedup_bam.out.bai
    deduplog = dedup_bam.out.logs
}


