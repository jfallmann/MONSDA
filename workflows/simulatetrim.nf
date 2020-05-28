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
    T1SAMPLES = T1 + T2
    T1SAMPLES.sort()

}else{
    T1SAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    T1SAMPLES.sort()
}

process simtrim{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("fastq.gz") > 0)               "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimmed.fastq.gz"
        else null
    }

    input:
    path read

    output:
    path "*.fastq.gz" , emit: trimmed_reads

    script:
    """
    echo $read
    """

}

workflow TRIMMING{
    take: dummy

    main:
    samples_ch1 = Channel.from(T1SAMPLES)

    simtrim(samples_ch1)

    emit:
    simtrim.out.trimmed_reads

}
