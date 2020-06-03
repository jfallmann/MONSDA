QCENV=params.QCENV ?: null
QCBIN=params.QCBIN ?: null

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
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $read
    """
}

workflow QC_MAPPING{
    take: dummy

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
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
        M1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+".bam"
        }
        M1SAMPLES.sort()
        mapped_samples_ch = Channel.fromPath(M1SAMPLES)

    }

    qc_mapped(mapped_samples_ch)

    emit:
    qc = qc_mapped.out.fastqc_results
}
