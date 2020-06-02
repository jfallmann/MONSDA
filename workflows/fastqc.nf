QCENV=params.QCENV ?: null
QCBIN=params.QCBIN ?: null

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
    samples_ch = Channel.fromPath(R1SAMPLES).merge(Channel.fromPath(R2SAMPLES))

    T1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1.fastq.gz"
    }
    T1SAMPLES.sort()
    T2SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2.fastq.gz"
    }
    T2SAMPLES.sort()
    trimmed_samples_ch = Channel.fromPath(T1SAMPLES).merge(Channel.fromPath(T2SAMPLES))

    M1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/"+element+"_R1_sorted.bam"
    }
    M1SAMPLES.sort()
    M2SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_sorted.bam"
    }
    M2SAMPLES.sort()
    mapped_samples_ch = Channel.fromPath(M1SAMPLES).merge(Channel.fromPath(M2SAMPLES))


}else{
    RSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    RSAMPLES.sort()
    samples_ch = Channel.fromPath(RSAMPLES)

        T1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+".fastq.gz"
    }
    T1SAMPLES.sort()
    trimmed_samples_ch = Channel.fromPath(T1SAMPLES)

    M1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/"+element+".bam"
    }
    M1SAMPLES.sort()
    mapped_samples_ch = Channel.fromPath(M1SAMPLES)

}



//PROCESSES
process qc_mapped{
    conda "${workflow.workDir}/../nextsnakes/envs/$QCENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    val dummy
    path read

    output:
    path "*.{zip,html}", emit: mapfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $read
    """
}

workflow QC_MAPPING{
    take: dummy

    main:

    qc_raw(dummy,samples_ch)
    qc_trimmed(dummy,trimmed_samples_ch)
    qc_mapped(dummy,mapped_samples_ch)

    emit:
    rawqc = qc_raw.out
    trimqc = qc_trimmed.out
    mapqc = qc_mapped.out
}
