QCENV=params.QCENV ?: null
QCBIN=params.QCBIN ?: null
QCPARAMS = params.fastqc_params_1 ?: ''

process collect_multi{
    input:
    path check

    output:
    path "collect.txt"

    script:
    """
    echo "STARTING MULTIQC" >> collect.txt
    """
}

process multiqc{
    conda "${workflow.workDir}/../NextSnakes/envs/$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1
    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$COMBO$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    path samples
    val collect

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -s .
    """
}

workflow MULTIQC{
    take: collection

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

        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        samples_ch.join(Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES)))

    }else{
        RSAMPLES=SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)

        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        samples_ch.join(Channel.fromPath(T1SAMPLES))
    }

    MSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted.bam"
    }
    MSAMPLES.sort()
    USAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+"_mapped_sorted_unique.bam"
    }
    USAMPLES.sort()

    samples_ch.join(Channel.fromPath(MSAMPLES))
    samples_ch.join(Channel.fromPath(USAMPLES))

    collect_multi(collection.collect())
    multiqc(samples_ch.collect(), collect_multi.out)

    emit:
    mqcres = multiqc.out.multiqc_results
}
