COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTIDX = get_always('COUNTINGIDX')
COUNTUIDX = get_always('COUNTINGUIDX')
COUNTUIDXNAME = get_always('COUNTINGUIDXNAME')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]
COUNTUIDX.replace('.idx','')

COUNTPARAMS = get_always('featurecounts_params_COUNT') ?: ''
FEAT = get_always('COUNTINGFEAT') ?: ''
COUNTMAP = get_always('COUNTINGMAP') ?: ''

//COUNTING PROCESSES
process count_fastq{
    conda "base.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"   
        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/countfastq.log"
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
        fn = file(reads[0]).getSimpleName().replaceAll(/\Q_R1\E/,"")    
        oo = fn+"_raw_R1_fq.count"        
        ft = file(reads[0]).getSimpleName().replaceAll(/\Q_R2\E/,"")    
        ot = ft+"_raw_R2_fq.count"
        ol = fn+".log"
        """
        a=$(zcat $r1|wc -l ); echo $((a/4)) > $oo;done 2>> {log} &&
        a=$(zcat $r2|wc -l ); echo $((a/4)) > $ot;done 2>> {log}
        """
    }else{
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q\E/,"")    
        oc = fn+"_raw_fq.count"
        ol = fn+".log"
        """
        a=$(zcat $read|wc -l ); echo $((a/4)) > $oc;done 2>> $ol
        """
    }
}

process count_trimmed_fastq{
    conda "base.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"   
        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/count_trimmedreads.log"
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
        fn = file(reads[0]).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")    
        oo = fn+"_raw_R1_fq.count"        
        ft = file(reads[0]).getSimpleName().replaceAll(/\Q_R2_trimmed\E/,"")    
        ot = ft+"_raw_R2_fq.count"
        ol = fn+".log"
        """
        a=$(zcat $r1|wc -l ); echo $((a/4)) > $oo;done 2>> {log} &&
        a=$(zcat $r2|wc -l ); echo $((a/4)) > $ot;done 2>> {log}
        """
    }else{
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")    
        oc = fn+"_raw_fq.count"
        ol = fn+".log"
        """
        a=$(zcat $read|wc -l ); echo $((a/4)) > $oc;done 2>> $ol
        """
    }
}

process count_mappers{
    conda "samtools.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"        
        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/count_mappers.log"

    }

    input:
    path reads

    output:
    path "*.count", emit: fq_cts

    script:        
    fn = file(reads).getSimpleName()
    oc = fn+"_mapped.count"
    ol = fn+".log"
    sortmem = '30%'
    """
    export LC_ALL=C; samtools view -F 260 $reads | cut -d$'\t' -f1|sort --parallel=$THREADS -S $sortmem -T TMP -u |wc -l > $oc;done 2>> $ol
    """
}

process featurecount{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "COUNTS/Featurecounts_$FEAT/$COMBO$CONDITION/${file(filename).getSimpleName()}.counts.gz"        
        else if (filename.indexOf(".summary") > 0)      "COUNTS/Featurecounts_$FEAT/$COMBO$CONDITION/${file(filename).getSimpleName()}.counts.summary"        
        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/featurecount_${FEAT}s.log"

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
    oc = fn+".counts.summary"
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
    $COUNTBIN -T $THREADS $COUNTPARAMS $pair $stranded $COUNTINGMAP -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=$THREADS -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc && mv tmpcts.summary $os
    """
}

process summarize_counts{
    conda "base.yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".summary") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.summary"   
        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/summarize_counts.log"
    }

    input:
    path reads

    output:
    path "*.summary", emit: sum
    path "*.log", emit: sum_log

    script:
    fn = file(reads).getSimpleName()
    oc = fn+".summary"
    ol = fn+".log"
    """
    echo -ne \"$reads\t\" >> $oc && if [[ -s $reads ]]; then cat $reads >> $oc; else echo '0' >> $oc 2> $ol
    """
}

workflow COUNTING{ //HIER WEITER
    take: collection

    main:
   
    checkidx = file(COUNTIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(COUNTUIDX)
        if (PAIRED == 'paired'){
            salmon_quant(idxfile.combine(samples_ch.collate(2)))
        } else{
            salmon_quant(idxfile.combine(samples_ch.collate(1)))
        }        
    }
    else{
        genomefile = Channel.fromPath(COUNTREF)
        salmon_idx(genomefile)
        if (PAIRED == 'paired'){
            salmon_quant(salmon_idx.out.idx.combine(samples_ch.collate(2)))
        } else{
            salmon_quant(salmon_idx.out.idx.combine(samples_ch.collate(1)))
        }
    }

    emit:
    counts = salmon_quant.out.counts
    logs = salmon_quant.out.logs
}

//process count_unique_mappers{
//    conda "samtools.yaml"
//    cpus THREADS
//    //validExitStatus 0,1
//
//    publishDir "${workflow.workDir}/../" , mode: 'link',
//    saveAs: {filename ->
//        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/count_unique_mappers.log"
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
//        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/count_dedup_mappers.log"
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
//        if (filename.indexOf(".count") > 0)      "COUNTS/$COMBO$CONDITION/${file(filename).getSimpleName()}.count"        
//        else if (filename.indexOf(".log") > 0)        "LOGS/$COMBO$CONDITION/${file(filename).getSimpleName()}/count_unique_dedup_mappers.log"
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