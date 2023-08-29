COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTIDX = get_always('COUNTINGIDX')
COUNTUIDX = get_always('COUNTINGUIDX')
COUNTUIDXNAME = get_always('COUNTINGUIDXNAME')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]
COUNTUIDX?.replace('.idx','') 

COUNTPARAMS = get_always('featurecounts_params_COUNT') ?: ''
FEAT = get_always('COUNTINGFEAT') ?: ''
COUNTMAP = get_always('COUNTINGMAP') ?: ''

//COUNTING PROCESSES
process count_fastq{
    conda "base.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"   
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/countfastq.log"
    }

    input:
    path reads

    output:
    path "*.count", emit: fq_cts
    path "*.log", emit: fq_log

    script:    
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1\E/,"")    
        oo = fn+"_raw_R1_fq.count"        
        ft = file(r2).getSimpleName().replaceAll(/\Q_R2\E/,"")    
        ot = ft+"_raw_R2_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $r1|wc -l ); echo \$((a/4)) > $oo 2>> $ol &&
        a=\$(zcat $r2|wc -l ); echo \$((a/4)) > $ot 2>> $ol
        """
    }else{
        fn = file(reads[0]).getSimpleName().replaceAll(/\Q\E/,"")    
        oc = fn+"_raw_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $reads|wc -l ); echo \$((a/4)) > $oc 2>> $ol
        """
    }
}

process count_trimmed_fastq{
    conda "base.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"   
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_trimmedreads.log"
    }

    input:
    path reads

    output:
    path "*.count", emit: fqt_cts
    path "*.log", emit: fqt_log

    script:    
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")    
        oo = fn+"_trimmed_R1_fq.count"        
        ft = file(r2).getSimpleName().replaceAll(/\Q_R2_trimmed\E/,"")    
        ot = ft+"_trimmed_R2_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $r1|wc -l ); echo \$((a/4)) > $oo 2>> $ol &&
        a=\$(zcat $r2|wc -l ); echo \$((a/4)) > $ot 2>> $ol
        """
    }else{
        fn = file(reads[0]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")    
        oc = fn+"_trimmed_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $reads|wc -l ); echo \$((a/4)) > $oc 2>> $ol
        """
    }
}

process count_dedup_fastq{
    conda "base.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"   
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_dedupreads.log"
    }

    input:
    path reads

    output:
    path "*.count", emit: fqd_cts
    path "*.log", emit: fqd_log

    script:    
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_dedup\E/,"")    
        oo = fn+"_dedup_R1_fq.count"        
        ft = file(r2).getSimpleName().replaceAll(/\Q_R2_dedup\E/,"")    
        ot = ft+"_dedup_R2_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $r1|wc -l ); echo \$((a/4)) > $oo 2>> $ol &&
        a=\$(zcat $r2|wc -l ); echo \$((a/4)) > $ot 2>> $ol
        """
    }else{
        fn = file(reads[0]).getSimpleName().replaceAll(/\Q_dedup\E/,"")    
        oc = fn+"_dedup_fq.count"
        ol = fn+".log"
        """
        a=\$(zcat $reads|wc -l ); echo \$((a/4)) > $oc 2>> $ol
        """
    }
}

process count_mappers{
    conda "samtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.count"        
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/count_mappers.log"

    }

    input:
    path reads

    output:
    path "*.count", emit: map_cts

    script:        
    fn = file(reads).getSimpleName()
    oc = fn+".count"
    ol = fn+".log"
    sortmem = '30%'
    """
    mkdir -p TMP; export LC_ALL=C; samtools view -F 260 $reads | cut -d\$'\\t' -f1|sort --parallel=${task.cpus} -S $sortmem -T TMP -u |wc -l > $oc 2>> $ol
    """
}

process featurecount{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/Featurecounts_${FEAT}s/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.counts.gz"        
        else if (filename.indexOf(".summary") > 0)      "COUNTS/Featurecounts_${FEAT}s/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}.counts.summary"        
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}/featurecount_${FEAT}s.log"

    }

    input:
    path fls

    output:
    path "*.counts.gz", emit: fc_cts
    path "*.summary", emit: fc_summary
    path "*.log", emit: fc_log

    script:
    anno = fls[0]
    reads = fls[1]        
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
    mkdir -p TMP; $COUNTBIN -T ${task.cpus} $COUNTPARAMS $pair $stranded $COUNTMAP -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=${task.cpus} -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc 2>> $ol && mv tmpcts.summary $os
    """
}

process summarize_counts{
    conda "base.yaml"
    cpus THREADS
	cache 'lenient'
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
    for i in $reads;do echo -ne \"\$i\\t\" >> summary && if [[ -s \$i ]]; then cat \$i >> summary; else echo '0' >> summary;fi; done 2>> log
    """
}

workflow COUNTING{ 
    take: collection

    main:

    if (PAIRED == 'paired'){
        RAWSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
        }
        TRIMSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"_{R2,R1}*.fastq.gz"
        }
        DEDUPSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/${COMBO}/"+element+"_{R2,R1}*.fastq.gz"
        }
    }else{
        RAWSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
        }
        TRIMSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"*.fastq.gz"
        }
        DEDUPSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/${COMBO}/"+element+"*.fastq.gz"
        }
    }
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"*.bam"
    }

    rawsamples_ch = Channel.fromPath(RAWSAMPLES.sort())
    trimsamples_ch = Channel.fromPath(TRIMSAMPLES.sort())
    dedupsamples_ch = Channel.fromPath(DEDUPSAMPLES.sort())
    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())  
    annofile = Channel.fromPath(COUNTANNO)
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
    featurecount(annofile.combine(mapsamples_ch.collate(1)))
    summarize_counts(count_fastq.out.fq_cts.concat(count_dedup_fastq.out.fqd_cts.concat(count_trimmed_fastq.out.fqt_cts.concat(count_mappers.out.map_cts))).collect())

    emit:
    counts = summarize_counts.out.sum
    logs = summarize_counts.out.sum_log
}