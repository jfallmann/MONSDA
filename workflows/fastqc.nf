QCENV=params.QCENV ?: null
QCBIN=params.QCBIN ?: null

MRSAMPLES = null

MRSAMPLES = LONGSAMPLES.collect{
    element -> return "${workflow.workDir}/../MAPPED/"+element+"_mapped_sorted.sam.gz"
}
MRSAMPLES.sort()

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
    samples_ch = Channel.from(FQSAMPLES)
    trsamples_ch = Channel.from(TRSAMPLES)
    mapsamples_ch = Channel.from(MRSAMPLES)

    qc_raw(samples_ch)
    qc_trimmed(trsamples_ch)
    qc_mapped(mapsamples_ch)

    emit:
    rawqc = qc_raw.out.fastqc_results
    trimqc = qc_trimmed.out.trfastqc_results
    mapqc = qc_mapped.out.mapfastqc_results
}
