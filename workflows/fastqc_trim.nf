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

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS --noextract -f fastq $read
    """
}

workflow QC_RAW{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_R{1,2}.*fastq.gz"
        }
        SAMPLES.sort()        
        samples_ch = Channel.fromPath(SAMPLES)//.join(Channel.fromPath(R2SAMPLES))

    }else{
        RSAMPLES=SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    collect_fqraw(collection.collect())
    qc_raw(collect_fqraw.out.done, samples_ch)
    //qc_raw(collection.collect())

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

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    //val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS --noextract -f fastq $read
    """
}

workflow QC_TRIMMING{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R{1,2}_trimmed.*fastq.gz"
        }
        SAMPLES.sort()        
        trimmed_samples_ch = Channel.fromPath(SAMPLES)//.join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    //collect_fqtrim(collection.collect())
    //qc_trimmed(collect_fqtrim.out.done, trimmed_samples_ch)
    qc_trimmed(collection.collect())

    emit:
    qc = qc_trimmed.out.fastqc_results
}
