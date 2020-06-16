//POST MAPPING PROCESSES

process collect_postmap{
    //echo true

    input:
    path dummy

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$dummy Collection successful!" > collect.txt
    """
}

process sortsam{
    conda "${workflow.workDir}/../nextsnakes/envs/samtools.yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".sam.gz") > 0)     "MAPPED/$CONDITION/$filename"
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
    set +o pipefail;samtools view -H $map|grep -P '^@HD' |pigz -p $THREADS -f > tmphead; samtools view -H $map|grep -P '^@SQ'|sort -t\$'\t' -k1,1 -k2,2V |pigz -p $THREADS -f >> tmphead ; samtools view -H $map|grep -P '^@RG'|pigz -p $THREADS -f >> tmphead ; samtools view -H $map|grep -P '^@PG'|pigz -p $THREADS -f >> tmphead ; export LC_ALL=C;zcat $map | grep -v \"^@\"|sort --parallel=$THREADS -S 25% -T TMP -t\$'\t' -k3,3V -k4,4n - |pigz -p $THREADS -f > tmpfile; cat tmphead tmpfile > $sorted && rm -f tmphead tmpfile $map && ln -s sortedsam.gz $map
    """
}

process sam2bam{
    conda "${workflow.workDir}/../nextsnakes/envs/samtools.yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".bam") > 0)     "MAPPED/$CONDITION/$filename"
        if (filename.indexOf(".log") > 0)     "MAPPED/$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path sam

    output:
    path "*_sorted.bam*", emit: bam
    path "*.log", emit: log

    script:
    fn = file(map[0]).getSimpleName()
    bam = fn+'_sorted.bam'
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads $THREADS && samtools index $bam 2> bam.log
    """
}

process uniqsam{
    conda "${workflow.workDir}/../nextsnakes/envs/samtools.yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("unique.sam.gz") > 0)     "UNIQUE_MAPPED/$CONDITION/$filename"
        if (filename.indexOf(".log") > 0)     "UNIQUE_MAPPED/$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path map

    output:
    path "*_sorted_unique.sam.gz", emit: sam
    path "*.log", emit: log

    script:
    fn = file(map[0]).getSimpleName()
    uniq = fn+'_sorted_unique.sam.gz'
    """
    $BINS/Shells/UniqueSam_woPicard.sh $map $uniq $THREADS 2> uniq.log
    """
}

process sam2bamuniq{
    conda "${workflow.workDir}/../nextsnakes/envs/samtools.yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".bam") > 0)     "UNIQUE_MAPPED/$CONDITION/$filename"
        if (filename.indexOf(".log") > 0)     "UNIQUE_MAPPED/$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path sam

    output:
    path "*_sorted_unique.bam*", emit: bam
    path "*.log", emit: log

    script:
    fn = file(map[0]).getSimpleName()
    bam = fn+'_sorted_unique.bam'
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads $THREADS && samtools index $bam 2> uniquebam.log
    """
}


workflow POSTMAPPING{
    take: samples_ch

    main:
    collect_postmap(samples_ch.collect())
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        M1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+"_R1_mapped.sam.gz"
        }
        M1SAMPLES.sort()
        M2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+"_R2_mapped.sam.gz"
        }
        M2SAMPLES.sort()
        mapped_samples_ch = Channel.fromPath(M1SAMPLES).merge(Channel.fromPath(M2SAMPLES))

    }else{
        M1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+"_mapped.sam.gz"
        }
        M1SAMPLES.sort()
        mapped_samples_ch = Channel.fromPath(M1SAMPLES)
    }


    sortsam(collect_postmap.out.done, mapped_samples_ch)
    sam2bam(collect_postmap.out.done, sortsam.out.sam)
    uniqsam(collect_postmap.out.done, sortsam.out.sam)
    sam2bamuniq(collect_postmap.out.done, uniqsam.out.sam)

    emit:
    postmap  = sam2bam.out.bam
    postmapuni = sam2bamuniq.out.bam
}
