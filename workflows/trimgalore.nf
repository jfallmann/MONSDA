TOOLENV=params.TRIMMINGENV ?: null
TOOLBIN=params.TRIMMINGBIN ?: null

TRIMPARAMS = params.trimgalore_params_0 ?: ''

T1SAMPLES = null
T2SAMPLES = null

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

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("{val_1.fq.gz,val_2.fq.gz}") > 0)               "TRIMMED_FASTQ/$CONDITION/${file(it).getSimpleName()}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)      "TRIMMED_FASTQ/$CONDITION/${file(it).getSimpleName()}trimming_report.txt"
        else null
    }

    input:
    path fread
    path sread

    output:
    path "*.fq.gz", emit: trimmed_reads
    path "*trimming_report.txt", emit: trim_report

    script:
    """
    $TOOLBIN --cores $THREADS --paired --gzip $TRIMPARAMS $fread $sread
    """
}

process trim{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("trimmed.fq.gz") > 0)              "TRIMMED_FASTQ/$CONDITION/${file(it).getSimpleName()}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)     "TRIMMED_FASTQ/$CONDITION/${file(it).getSimpleName()}trimming_report.txt"
        else null
    }

    input:
    path read

    output:
    path "*trimmed.fastq.gz", emit: trimmed_reads
    path "*trimming_report.txt", emit: trim_report

    script:
    """
    $TOOLBIN --cores $THREADS --gzip $TRIMPARAMS $read
    """
}

workflow TRIMMING{
    take: dummy

    main:
    samples_ch1 = Channel.from(T1SAMPLES)
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
