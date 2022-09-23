BINS = get_always('BINS')
DEENV = get_always('DEENV')
DEBIN = get_always('DEBIN')
DEREF = get_always('DEREF')
DEREFDIR = get_always('DEREFDIR')
DEANNO = get_always('DEANNO')
COUNTPARAMS = get_always('deseq2_params_COUNT') ?: ''
DEPARAMS = get_always('deseq2_params_DE') ?: ''
DEREPS = get_always('DEREPS') ?: ''
DECOMP = get_always('DECOMP') ?: ''
DECOMPS = get_always('DECOMPS') ?: ''
PCOMBO = get_always('COMBO') ?: 'none'

COUNTBIN, COUNTENV = ['featureCounts','countreads_de']
//DE PROCESSES

process featurecount{
    conda "$DEENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "DE/${SCOMBO}Featurecounts/${file(filename).getSimpleName()}.counts.gz"                
        else if (filename.indexOf(".log") > 0)        "LOGS/DE/${SCOMBO}/${file(filename).getSimpleName()}/featurecounts_deseq2_unique.log"
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
    $COUNTBIN -T $THREADS $COUNTPARAMS $pair $stranded -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=$THREADS -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc 2>> $ol && mv tmpcts.summary $os
    """
}

process prepare_count_table{
    conda "$DEENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "${COMBO}_COUNTS.gz")      "DE/${SCOMBO}/Tables/${COMBO}_COUNTS.gz"
        else if (filename == "${COMBO}_ANNOTATION.gz")      "DE/${SCOMBO}/Tables/${COMBO}_ANNOTATION.gz"
        else if (filename == "log")      "LOGS/DE/${SCOMBO}/${COMBO}_prepare_count_table.log"
    }

    input:
    //path '*.count*'// from reads
    path reps

    output: 
    path "*_COUNTS.gz", emit: counts
    path "*_ANNOTATION.gz", emit: anno
    path "log", emit: log

    script:
    """
    ${BINS}/Analysis/build_count_table.py $DEREPS --table ${COMBO}_COUNTS.gz --anno ${COMBO}_ANNOTATION.gz 2> log
    """
}

process run_deseq2{
    conda "$DEENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DE/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("_figure") > 0)      "DE/${SCOMBO}/Figures/${file(filename).getName()}"                
        else if (filename.indexOf("SESSION") > 0)      "DE/${SCOMBO}/${file(filename).getName()}"                
        else if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${SCOMBO}.Rmd"                
        else if (filename.indexOf("log") > 0)        "LOGS/DE/${SCOMBO}/run_deseq2.log"
    }

    input:
    //path '*.count*'// from reads
    path cts
    path anno

    output:
    path "summary", emit: sum
    path "log", emit: sum_log

    script:    
    outdir = "DE"+File.separatorChar+${SCOMBO}
    """
    Rscript --no-environ --no-restore --no-save $BINS $anno $cts $DEANNO $outdir $DECOMP $PCOMBO $THREADS $DEPARAMS 2> log
    """
}

process filter_significant{
    
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