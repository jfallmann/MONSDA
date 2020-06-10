MAPENV=params.MAPPINGENV ?: null
MAPBIN=params.MAPPINGBIN ?: null

IDXPARAMS = params.star_params_0 ?: ''
MAPPARAMS = params.star_params_1 ?: ''

//PROCESSES
process create_idx{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".fq.gz") > 0)      "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)     "TRIMMED_FASTQ/$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else null
    }

    input:
    //val collect
    path reads

    output:
    path "*{val,trimmed}*.fq.gz", emit: trim
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

workflow MAPPING{
    take: samples_ch

    main:
    //collect_results(samples_ch.collect())
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).merge(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    star_idx(samples_ch)
    star(star_idx.out.idx, samples_ch)

    emit:
    idx = star_idx.out.idx
    mapped  = star.out.maps
}
