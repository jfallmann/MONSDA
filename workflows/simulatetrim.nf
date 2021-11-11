T1SAMPLES = null
T2SAMPLES = null

process simtrim{
    //conda "$TOOLENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_trimmed.fastq.gz") > 0)           "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)     "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    path read

    output:
    path "*trimmed.fastq.gz" , emit: trim
    path "*trimming_report.txt", emit: rep

    script:
    a=read.getSimpleName()+"_trimming_report.txt"
    b=read.getSimpleName()+"_trimmed.fastq.gz"
    """
    ln -sf $read $b ; echo "simulated $read trimming" > $a
    """

}

workflow TRIMMING{
    take: samples_ch

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        if (RUNDEDUP == 'enabled'){
            R1SAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_R1_dedup.fastq.gz"
            }
            R1SAMPLES.sort()
            R2SAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_R2_dedup.fastq.gz"
            }
            R2SAMPLES.sort()            
        }
        else{   
            R1SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
            }
            R1SAMPLES.sort()
            R2SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
            }
            R2SAMPLES.sort()            
        }
        samples_ch = Channel.fromPath(R1SAMPLES).join(Channel.fromPath(R2SAMPLES))
    }else{
        if (RUNDEDUP == 'enabled'){
            RSAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_dedup.fastq.gz"
            }
        }
        else{
            RSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
            }
        }                 
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    simtrim(samples_ch)

    emit:
    trimmed = simtrim.out.trim
    report  = simtrim.out.rep
}