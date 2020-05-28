TRIMENV=params.TRIMMINGENV ?: null
TRIMBIN=params.TRIMMINGBIN ?: null

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
    conda "${workflow.workDir}/../nextsnakes/envs/$TRIMENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("fq.gz") > 0)               "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)      "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    path fread
    path sread

    output:
    path "*val_1.fq.gz", emit: trimmed_r1
    path "*val_2.fq.gz", emit: trimmed_r2
    path "*trimming_report.txt"

    script:
    """
    $TRIMBIN --cores $THREADS --paired --gzip $TRIMPARAMS $fread $sread
    """
}

process trim{
    conda "${workflow.workDir}/../nextsnakes/envs/$TRIMENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("trimmed.fq.gz") > 0)      "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf("report.txt") >0)     "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    path read

    output:
    path "*trimmed.fq.gz", emit: trimmed
    path "*trimming_report.txt"

    script:
    """
    $TRIMBIN --cores $THREADS --gzip $TRIMPARAMS $read
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
        trimmed_r1 = trim_paired.out.trimmed_r1
        trimmed_r2 = trim_paired.out.trimmed_r2
    }
    else{
        trim(samples_ch1)

        emit:
        trimmed = trim.out.trimmed
    }
}
