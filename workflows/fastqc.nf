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
    val collect
    path map
    path uni

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f sam_mapped $map $uni
    """
}

workflow QC_MAPPING{
    take: mapped_samples_ch

    main:
    x = collect_results(mapped_samples_ch.collect())
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        M1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+"_R1_mapped_sorted.bam"
        }
        M1SAMPLES.sort()
        M2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+"_R2_mapped_sorted.bam"
        }
        M2SAMPLES.sort()
        U1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../UNIQUE_MAPPED/"+element+"_R1_mapped_sorted_unique.bam"
        }
        U1SAMPLES.sort()
        U2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../UNIQUE_MAPPED/"+element+"_R2_mapped_sorted_unique.bam"
        }
        U2SAMPLES.sort()

        mapped_samples_ch = Channel.fromPath(M1SAMPLES).merge(Channel.fromPath(M2SAMPLES))
        unique_samples_ch = Channel.fromPath(U1SAMPLES).merge(Channel.fromPath(U2SAMPLES))


    }else{
        M1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/"+element+".bam"
        }
        M1SAMPLES.sort()
        U1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../UNIQUE_MAPPED/"+element+".bam"
        }
        U1SAMPLES.sort()

        mapped_samples_ch = Channel.fromPath(M1SAMPLES)
        unique_samples_ch = Channel.fromPath(U1SAMPLES)

    }

    qc_mapped(x.out.done, mapped_samples_ch, unique_samples_ch)

    emit:
    qc = qc_mapped.out.fastqc_results
}
