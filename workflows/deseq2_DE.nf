DEENV = get_always('DEENV')
DEBIN = get_always('DEBIN')
DEREF = get_always('DEREF')
DEREFDIR = get_always('DEREFDIR')
DEANNO = get_always('DEANNO')

DEPARAMS = get_always('featurecounts_params_COUNT') ?: ''
FEAT = get_always('DEFEAT') ?: ''
DEMAP = get_always('DEMAP') ?: ''

//DE PROCESSES

process featurecount{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/Featurecounts_$FEAT/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.counts.gz"        
        else if (filename.indexOf(".summary") > 0)      "COUNTS/Featurecounts_$FEAT/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.counts.summary"        
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/featurecount_${FEAT}s.log"

    }

    input:
    path anno
    path reads

    output:
    path "*.counts.gz", emit: fc_cts
    path "*.summary", emit: fc_summary
    path "*.log", emit: fc_log

    script:        
    fn = file(reads).getSimpleName()
    oc = fn+".counts.gz"
    os = fn+".counts.summary"
    ol = fn+".log"
    sortmem = '30%'
    if (PAIRED == 'paired'){
        pair = "-p"
    }
    else{
        pair= ""
    }
    if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '-s 1'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
            stranded = '-s 2'
        }else{
            stranded = ''
    }
    """
    $COUNTBIN -T $THREADS $COUNTPARAMS $pair $stranded $COUNTMAP -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=$THREADS -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc 2>> $ol && mv tmpcts.summary $os
    """
}

process summarize_counts{
    conda "base.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "summary")      "COUNTS/${SCOMBO}/${CONDITION}/summary"
        else if (filename == "log")        "LOGS/${SCOMBO}/${CONDITION}/summarize_counts.log"
    }

    input:
    //path '*.count*'// from reads
    path reads

    output:
    path "summary", emit: sum
    path "log", emit: sum_log

    script:    
    """
    for i in $reads;do echo -ne \"\$i\t\" >> summary && if [[ -s \$i ]]; then cat \$i >> summary; else echo '0' >> summary;fi; done 2>> log
    """
}

workflow DE{ 
    take: collection

    main:

    if (PAIRED == 'paired'){
        RAWSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
        }
        TRIMSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}"+element+"_{R2,R1}*.fastq.gz"
        }
        DEDUPSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/${COMBO}"+element+"_{R2,R1}*.fastq.gz"
        }
    }else{
        RAWSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
        }
        TRIMSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}"+element+"*.fastq.gz"
        }
        DEDUPSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/${COMBO}"+element+"*.fastq.gz"
        }
    }
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}"+element+"*.bam"
    }

    rawsamples_ch = Channel.fromPath(RAWSAMPLES)
    trimsamples_ch = Channel.fromPath(TRIMSAMPLES)
    dedupsamples_ch = Channel.fromPath(DEDUPSAMPLES)
    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES)  
    annofile = Channel.fromPath(DEANNO)
    //annofile.subscribe {  println "ANNO: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }


    if (PAIRED == 'paired'){
        count_fastq(rawsamples_ch.collate(2))
        count_trimmed_fastq(trimsamples_ch.collate(2))
        count_dedup_fastq(dedupsamples_ch.collate(2))
    } else{
        count_fastq(rawsamples_ch.collate(1))
        count_trimmed_fastq(trimsamples_ch.collate(1))
        count_dedup_fastq(dedupsamples_ch.collate(1))
    }        
    count_mappers(mapsamples_ch.collate(1))
    featurecount(annofile, mapsamples_ch.collate(1))
    summarize_counts(count_fastq.out.fq_cts.concat(count_dedup_fastq.out.fqd_cts.concat(count_trimmed_fastq.out.fqt_cts.concat(count_mappers.out.map_cts))).collect())

    emit:
    counts = summarize_counts.out.sum
    logs = summarize_counts.out.sum_log
}

//process count_unique_mappers{
//    conda "samtools.yaml"
//    cpus THREADS
//    //validExitStatus 0,1
//
//    publishDir "${workflow.workDir}/../" , mode: 'link',
//    saveAs: {filename ->
//        if (filename.indexOf(".count") > 0)      "COUNTS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_unique_mappers.log"
//
//    }
//
//    input:
//    path reads
//
//    output:
//    path "*.count", emit: fq_cts
//
//    script:        
//    fn = file(reads).getSimpleName()
//    oc = fn+"_mapped_unique.count"
//    ol = fn+".log"
//    sortmem = '30%'
//    """
//    export LC_ALL=C; samtools view -F 260 $reads | cut -d$'\t' -f1|sort --parallel=$THREADS -S $sortmem -T TMP -u |wc -l > $oc;done 2>> $ol
//    """
//}
//
//process count_dedup_mappers{
//    conda "samtools.yaml"
//    cpus THREADS
//    //validExitStatus 0,1
//
//    publishDir "${workflow.workDir}/../" , mode: 'link',
//    saveAs: {filename ->
//        if (filename.indexOf(".count") > 0)      "COUNTS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_dedup_mappers.log"
//
//    }
//
//    input:
//    path reads
//
//    output:
//    path "*.count", emit: fq_cts
//
//    script:        
//    fn = file(reads).getSimpleName()
//    oc = fn+"_mapped_dedup.count"
//    ol = fn+".log"
//    sortmem = '30%'
//    """
//    export LC_ALL=C; samtools view -F 260 $reads | cut -d$'\t' -f1|sort --parallel=$THREADS -S $sortmem -T TMP -u |wc -l > $oc;done 2>> $ol
//    """
//}
//
//process count_unique_dedup_mappers{
//    conda "samtools.yaml"
//    cpus THREADS
//    //validExitStatus 0,1
//
//    publishDir "${workflow.workDir}/../" , mode: 'link',
//    saveAs: {filename ->
//        if (filename.indexOf(".count") > 0)      "COUNTS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_unique_dedup_mappers.log"
//
//    }
//
//    input:
//    path reads
//
//    output:
//    path "*.count", emit: fq_cts
//
//    script:        
//    fn = file(reads).getSimpleName()
//    oc = fn+"_mapped_unique_dedup.count"
//    ol = fn+".log"
//    sortmem = '30%'
//    """
//    export LC_ALL=C; samtools view -F 260 $reads | cut -d$'\t' -f1|sort --parallel=$THREADS -S $sortmem -T TMP -u |wc -l > $oc;done 2>> $ol
//    """
//}