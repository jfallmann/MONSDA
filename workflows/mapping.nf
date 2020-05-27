if (PAIRED == 'paired'){
    M1 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/TRIMMED_FASTQ/"+element+"_R1.fastq.gz"
    }
    M2 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/TRIMMED_FASTQ/"+element+"_R2.fastq.gz"
    }
    M1SAMPLES = M1.sort()
    M2SAMPLES = M2.sort()

}else{
    M1SAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/FASTQ/"+element+".fastq.gz"
    }
    M1SAMPLES.sort()
}
