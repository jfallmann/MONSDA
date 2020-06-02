QCRENV=params.QCENV ?: null
QCRBIN=params.QCBIN ?: null


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

}else{
    RSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    RSAMPLES.sort()
    samples_ch = Channel.fromPath(RSAMPLES)
}

process qc_raw{
    conda "${workflow.workDir}/../nextsnakes/envs/$QCRENV"+".yaml"
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

workflow QC_RAW{
    take: dummy

    main:
    qc_raw(samples_ch)

    emit:
    rawqc = qc_raw.out
}
