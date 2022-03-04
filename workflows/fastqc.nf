QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_QC') ?: ''

//QC RAW
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
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS --noextract -f fastq $read
    """
}

workflow QC_RAW{
    take:
    collection

    main:
    
    qc_raw(samples_ch)

    emit:
    qc = qc_raw.out.fastqc_results
}

//QC TRIM

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
    
    qc_trimmed(collection)

    emit:
    qc = qc_trimmed.out.fastqc_results
}


//QC MAP

process qc_mapped{
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
    path map

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS -f bam $map
    """
}

workflow QC_MAPPING{
    take: collection
    main:
   
    qc_mapped(collection)

    emit:
    qc = qc_mapped.out.fastqc_results
}

// DEDUP QC

process qc_dedup{
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
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS --noextract -f fastq $read
    """
}

workflow QC_DEDUP{
    take: collection

    main:
   
    qc_dedup(collection)

    emit:
    qc = qc_dedup.out.fastqc_results
}
