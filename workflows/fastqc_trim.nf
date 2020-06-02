QCTENV=params.QCENV ?: null
QCTBIN=params.QCBIN ?: null

//SAMPLE CHANNELS
if (PAIRED == 'paired'){
    T1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1.fastq.gz"
    }
    T1SAMPLES.sort()
    T2SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2.fastq.gz"
    }
    T2SAMPLES.sort()
    trimmed_samples_ch = Channel.fromPath(T1SAMPLES).merge(Channel.fromPath(T2SAMPLES))

}else{
    T1SAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+".fastq.gz"
    }
    T1SAMPLES.sort()
    trimmed_samples_ch = Channel.fromPath(T1SAMPLES)

}


process qc_trimmed{
    conda "${workflow.workDir}/../nextsnakes/envs/$QCTENV"+".yaml"
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

workflow QC_TRIMMING{
    take: dummy

    main:
    qc_trimmed(trimmed_samples_ch)

    emit:
    trimqc = qc_trimmed.out
}
