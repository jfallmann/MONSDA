MRSAMPLES = null

MRSAMPLES = LONGSAMPLES.collect{
    element -> return "${workflow.workDir}/../MAPPED/"+element+"_mapped_sorted.sam.gz"
}
MRSAMPLES.sort()

process qc_mapped{
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
    path mapped

    output:
    path "*.{zip,html}", emit: mapfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $read
    """
}

workflow QC_MAPPING{
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
