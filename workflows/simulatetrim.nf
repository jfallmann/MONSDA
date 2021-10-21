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
    conda "$TOOLENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("fastq.gz") > 0)           "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)     "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    path read

    output:
    path "*trimmed.fastq.gz" , emit: trim
    path "*trimming_report.txt", emit: rep

    script:
    """
    a=basename($read); rename 's/.fastq.gz/trimmed.fastq.gz/' $read ; echo "simulated" > \${a}_trimming_report.txt
    """

}

workflow TRIMMING{
    take: samples_ch

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        R1SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
        }
        R1SAMPLES.sort()
        R2SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
        }
        R2SAMPLES.sort()
        samples_ch = Channel.fromPath(R1SAMPLES).join(Channel.fromPath(R2SAMPLES))
    }else{
        RSAMPLES=SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    simtrim(samples_ch)

    emit:
    trimmed = simtrim.out.trim
    report  = simtrim.out.rep
}
