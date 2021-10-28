QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_QC') ?: ''

//QC RAW
process collect_fqraw{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process qc_raw{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

workflow QC_RAW{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        R1SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
        }
        R1SAMPLES.sort()
        R2SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
        }
        R2SAMPLES.sort()
        samples_ch = Channel.fromPath(R1SAMPLES).join(Channel.fromPath(R2SAMPLES))

    }else{
        RSAMPLES=SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    collect_fqraw(collection.collect())
    qc_raw(collect_fqraw.out.done, samples_ch)

    emit:
    qc = qc_raw.out.fastqc_results
}

//QC TRIM

process collect_fqtrim{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process qc_trimmed{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

workflow QC_TRIMMING{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    collect_fqtrim(collection.collect())
    qc_trimmed(collect_fqtrim.out.done, trimmed_samples_ch)

    emit:
    qc = qc_trimmed.out.fastqc_results
}


//QC MAP

process collect_fqmap{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}


//PROCESSES
process qc_mapped{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path map
    path uni

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f bam_mapped $map $uni
    """
}

workflow QC_MAPPING{
    take: collection

    main:
    //SAMPLE CHANNELS
    M1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted.bam"
    }
    M1SAMPLES.sort()
    U1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted_unique.bam"
    }
    U1SAMPLES.sort()

    mapped_samples_ch = Channel.fromPath(M1SAMPLES)
    unique_samples_ch = Channel.fromPath(U1SAMPLES)

    collect_fqmap(collection.collect())
    qc_mapped(collect_fqmap.out.done, mapped_samples_ch, unique_samples_ch)

    emit:
    qc = qc_mapped.out.fastqc_results
}

// DEDUP QC

process collect_fqdedup{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process qc_dedup{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

workflow QC_DEDUP{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_R1_dedup.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_R2_dedup.fastq.gz"
        }
        T2SAMPLES.sort()
        dedup_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_dedup.fastq.gz"
        }
        T1SAMPLES.sort()
        dedup_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    collect_fqdedup(collection.collect())
    qc_dedup(collect_fqdedup.out.done, dedup_samples_ch)

    emit:
    qc = qc_dedup.out.fastqc_results
}
