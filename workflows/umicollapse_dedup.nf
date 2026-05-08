DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

DEDUPPARAMS = get_always('umicollapse_params_DEDUP') ?: ''

process dedup_bam{
    conda "$DEDUPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$DEDUPENV"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))          "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_dedup.bam.bai") > 0) "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("dedup.log") > 0)           "LOGS/${COMBO}/${CONDITION}/DEDUP/${file(filename).getName()}"
        else null
    }

    input:
    path todedup
    path bami
        
    output:
    path "*_dedup.bam", emit: bam
    path "*_dedup.bam.bai", emit: bai
    path "*_dedup.log", emit: logs

    memory { 20.GB * (1 << ((task.attempt ?: 1) - 1)) }
    time { 8.h * (1 << ((task.attempt ?: 1) - 1)) }

    script:
    bams = todedup[0]
    bais = todedup[1]
    outf = bams.getSimpleName()+"_dedup.bam"
    outl = bams.getSimpleName()+"_dedup.log"
    if (PAIRED == 'paired'){        
        """
            mkdir -p TMP && $DEDUPBIN bam -i $bams -o TMP/dedup.unsorted.bam --paired $DEDUPPARAMS > $outl 2>&1 && samtools sort -@ ${task.cpus} -o $outf TMP/dedup.unsorted.bam >> $outl 2>&1 && samtools index $outf >> $outl 2>&1
        """
    }
    else{
        """
            mkdir -p TMP && $DEDUPBIN bam -i $bams -o TMP/dedup.unsorted.bam $DEDUPPARAMS > $outl 2>&1 && samtools sort -@ ${task.cpus} -o $outf TMP/dedup.unsorted.bam >> $outl 2>&1 && samtools index $outf >> $outl 2>&1
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
