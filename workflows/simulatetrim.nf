T1SAMPLES = null
T2SAMPLES = null

process trim{
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
    take: 
    collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        if (RUNDEDUP == 'enabled'){
            SAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_R{1,2}_dedup.*fastq.gz"
            }
            SAMPLES.sort()            
        }
        else{   
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_R{1,2}.*fastq.gz"
            }
            SAMPLES.sort()            
        }
        samples_ch = Channel.fromPath(SAMPLES)//.join(Channel.fromPath(R2SAMPLES))
    }else{
        if (RUNDEDUP == 'enabled'){
            SAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_dedup.fastq.gz"
            }
        }
        else{
            SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
            }
        }                 
        SAMPLES.sort()
        samples_ch = Channel.fromPath(SAMPLES)
    }
    //samples_ch.mix(collection.collect())  //NEED FIX HERE, EITHER EMPTY OR NOT ABSPATH
    //simtrim(samples_ch.collect())
    collection.collect().join(samples_ch).unique().filter(!~/MONSDA.log/)
    trim(collection.collect())

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}