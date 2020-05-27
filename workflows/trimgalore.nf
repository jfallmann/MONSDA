TOOLENV=params.TRIMMINGENV ?: null
TOOLBIN=params.TRIMMINGBIN ?: null

TRIMPARAMS = params.trimgalore_params_0 ?: ''

if (PAIRED == 'paired'){
    T1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    T2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    T1SAMPLES = T1.sort()
    T2SAMPLES = T2.sort()

}else{
    T1SAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    T1SAMPLES.sort()
}

process trim_paired{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'move',
    saveAs: {filename ->
        if (filename.indexOf("val_*fq.gz") > 0)          "TRIMMED_FASTQ/$CONDITION/${filename.baseName}"+"_trimmed.fastq.gz"
        else null
    }

    input:
    path fread
    path sread

    output:
    path "*_trimmed.fastq.gz", emit: trimmed_reads

    script:
    """
    $TOOLBIN --cores $THREADS --paired --no_report_file --gzip $TRIMPARAMS $fread $sread
    """
}

process trim{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'move',
    saveAs: {filename ->
        if (filename.indexOf("*fq.gz") > 0)          "TRIMMED_FASTQ/$CONDITION/${filename.baseName}"+"_trimmed.fastq.gz"
        else null
    }

    input:
    path read

    output:
    path "*_trimmed.fastq.gz", emit: trimmed_reads

    script:
    """
    $TOOLBIN --cores $THREADS --no_report_file --gzip $TRIMPARAMS $read
    """
}

workflow TRIMMING{
    samples_ch1 = Channel.from(T1SAMPLES)
    take: bla ?: null
    main:
    if (PAIRED == 'paired'){
        samples_ch2 = Channel.from(T2SAMPLES)
        trim_paired(samples_ch1, samples_ch2)
        emit:
        trim_paired.out.trimmed_reads

    }
    else{
        trim(samples_ch1)
        emit:
        trim.out.trimmed_reads
    }
}
