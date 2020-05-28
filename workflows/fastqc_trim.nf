QCTENV=params.QCENV ?: null
QCTBIN=params.QCBIN ?: null

TRSAMPLES = null

if (PAIRED == 'paired'){
    TR1 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
    }
    TR2 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
    }
    TRSAMPLES = TR1+TR2
    TRSAMPLES.sort()

}else{
    TRSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
    }
    TRSAMPLES.sort()
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
    trsamples_ch = Channel.from(TRSAMPLES)
    qc_trimmed(trsamples_ch)

    emit:
    trimqc = qc_trimmed.out.fastqc_results

}
