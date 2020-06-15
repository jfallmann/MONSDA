TRIMENV=params.TRIMMINGENV ?: null
TRIMBIN=params.TRIMMINGBIN ?: null

TRIMPARAMS = params.trimgalore_params_0 ?: ''

//PROCESSES
process trim{
    conda "${workflow.workDir}/../nextsnakes/envs/$TRIMENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("R1_val_1.fq.gz") > 0)                  "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("R2_val_2.fq.gz") > 0)             "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("_trimmed.fq.gz") > 0)             "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)                  "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    //val collect
    path reads

    output:
    path "*.fq.gz", emit: trim
    path "*trimming_report.txt", emit: rep

    script:
    if (PAIRED == 'paired'){
        """
        $TRIMBIN --cores $THREADS --paired --gzip $TRIMPARAMS $reads
        """
    }
    else{
        """
        $TRIMBIN --cores $THREADS --gzip $TRIMPARAMS $reads
        """
    }
}

workflow TRIMMING{
    take: samples_ch

    main:
    //collect_results(samples_ch.collect())
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
        samples_ch = Channel.fromPath(R1SAMPLES).merge(Channel.fromPath(R2SAMPLES))
    }else{
        RSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    trim(samples_ch)

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}
