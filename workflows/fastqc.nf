FQSAMPLES = null
TRSAMPLES = null
MRSAMPLES = null


if (PAIRED == 'paired'){
    FR1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    FR2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    FQSAMPLES = FR1+FR2
    FQSAMPLES.sort()

    TR1 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
    }
    TR2 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
    }
    TRSAMPLES = TR1+TR2
    TRSAMPLES.sort()

}else{
    FQSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    FQSAMPLES.sort()

    TRSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
    }
    TRSAMPLES.sort()
}

MRSAMPLES = LONGSAMPLES.collect{
    element -> return "${workflow.workDir}/../MAPPED/"+element+"_mapped_sorted.sam.gz"
}
MRSAMPLES.sort()

process qc_raw{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

process qc_trimmed{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path read

    output:
    path "*.{zip,html}", emit: trfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

process qc_mapped{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path mapped

    output:
    path "*.{zip,html}", emit: mapfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $read
    """
}

workflow fastqc{
    samples_ch = Channel.from(FQSAMPLES)
    trsamples_ch = Channel.from(TRSAMPLES)
    mapsamples_ch = Channel.from(MRSAMPLES)

    main:
    qc_raw(samples_ch)
    qc_trimmed(trsamples_ch)
    qc_mapped(mapsamples_ch)

    emit:
    qc_raw.out.fastqc_results
    qc_trimmed.out.trfastqc_results
    qc_mapped.out.mapfastqc_results
}
