//POST MAPPING PROCESSES

process collect_postmap{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process sortsam{
    conda "samtools.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".sam.gz") > 0)     "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.sam.gz"
        else null
    }

    input:
    val collect
    path map

    output:
    path "*_sorted.sam.gz", emit: sam

    script:
    fn = file(map[0]).getSimpleName()
    sorted = fn+'_sorted.sam.gz'
    """
    set +o pipefail;samtools view -H $map|grep -P '^@HD' |pigz -p $THREADS -f > tmphead; samtools view -H $map|grep -P '^@SQ'|sort -t\$'\t' -k1,1 -k2,2V |pigz -p $THREADS -f >> tmphead ; samtools view -H $map|grep -P '^@RG'|pigz -p $THREADS -f >> tmphead ; samtools view -H $map|grep -P '^@PG'|pigz -p $THREADS -f >> tmphead ; export LC_ALL=C;zcat $map | grep -v \"^@\"|sort --parallel=$THREADS -S 25% -T TMP -t\$'\t' -k3,3V -k4,4n - |pigz -p $THREADS -f > tmpfile; cat tmphead tmpfile > $sorted && rm -f tmphead tmpfile ${workflow.workDir}/../MAPPED/$COMBO$CONDITION/$map
    """
}

process sam2bam{
    conda "samtools.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".bam") > 0)       "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bai"
        else if (filename.indexOf(".log") > 0)  "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    val collect
    path sam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: log

    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+".bam"
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads $THREADS && samtools index $bam 2> sam2bam.log
    """
}

process uniqsam{
    conda "samtools.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("unique.sam.gz") > 0)   "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.sam.gz"
        else if (filename.indexOf(".log") > 0)       "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    val collect
    path sam

    output:
    path "*_unique.sam.gz", emit: sam
    path "*.log", emit: log

    script:
    fn = file(sam[0]).getSimpleName()
    uniq = fn+'_unique.sam.gz'
    """
    $BINS/Shells/UniqueSam_woPicard.sh $sam $uniq $THREADS 2> uniq.log
    """
}

process sam2bamuniq{
    conda "samtools.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".bam") > 0)       "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.bai"
        else if (filename.indexOf(".log") > 0)  "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    val collect
    path sam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: log

    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+'.bam'
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads $THREADS && samtools index $bam 2> uniquebam.log
    """
}


workflow POSTMAPPING{
    take: collection

    main:
    //SAMPLE CHANNELS
    M1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped.sam.gz"
    }
    M1SAMPLES.sort()
    mapped_samples_ch = Channel.fromPath(M1SAMPLES)

    collect_postmap(collection.collect())
    sortsam(collect_postmap.out.done, mapped_samples_ch)
    sam2bam(collect_postmap.out.done, sortsam.out.sam)
    uniqsam(collect_postmap.out.done, sortsam.out.sam)
    sam2bamuniq(collect_postmap.out.done, uniqsam.out.sam)

    emit:
    postmap  = sam2bam.out.bam
    postmapuni = sam2bamuniq.out.bam
}
