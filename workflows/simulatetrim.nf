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
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_{R1,R2}_dedup.*fastq.gz"
            }    
        }
        else{   
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R1,R2}.*fastq.gz"
            } 
        }
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
    }
   if (collection.collect().contains('MONSDA.log')){
        if (PAIRED == 'paired'){
            collection = Channel.fromFilePairs(SAMPLES)
        }
        else{
            collection = Channel.fromPath(SAMPLES)
        }
    }
    trim(collection.collect())

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}