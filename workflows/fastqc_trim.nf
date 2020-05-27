TOOLENV=params.QCENV ?: null
TOOLBIN=params.QCBIN ?: null

FQSAMPLES = null
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
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
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
    path "*.{zip,html}", emit: trfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

workflow QC_TRIMMING{
    samples_ch = Channel.from(FQSAMPLES)
    trsamples_ch = Channel.from(TRSAMPLES)
    take: bla ?: null
    main:
    qc_trimmed(trsamples_ch)

    emit:
    qc_trimmed.out.trfastqc_results
}
