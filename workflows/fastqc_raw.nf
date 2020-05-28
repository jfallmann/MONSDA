QCRENV=params.QCENV ?: null
QCRBIN=params.QCBIN ?: null

FQSAMPLES = null

if (PAIRED == 'paired'){
    R1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    R2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    FQSAMPLES = R1+R2
    FQSAMPLES.sort()

}else{
    FQSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    FQSAMPLES.sort()
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
    samples_ch = Channel.from(FQSAMPLES)
    qc_raw(samples_ch)

    emit:
    rawqc = qc_raw.out
}
