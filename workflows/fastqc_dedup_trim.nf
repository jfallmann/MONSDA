QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_QC') ?: ''

// RAW QC
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

// TRIMMED QC

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
