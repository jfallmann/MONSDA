QCTENV=params.QCENV ?: null
QCTBIN=params.QCBIN ?: null

QCPARAMS = params.fastqc_params_0 ?: ''

process collect_fqtrim{
    //echo true

    input:
    path dummy

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$dummy Collection successful!" > collect.txt
    """
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
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

workflow QC_TRIMMING{
    take: trimmed_samples_ch

    main:
    collect_fqtrim(trimmed_samples_ch.collect())
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).merge(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    qc_trimmed(collect_fqtrim.out.done, trimmed_samples_ch)

    emit:
    qc = qc_trimmed.out.fastqc_results
}
