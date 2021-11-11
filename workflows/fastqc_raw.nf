QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_QC') ?: ''

process collect_fqraw{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process qc_raw{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    val collect
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS $QCPARAMS --noextract -f fastq $read
    """
}

workflow QC_RAW{
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

    }else{
        RSAMPLES=SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    collect_fqraw(collection.collect())
    qc_raw(collect_fqraw.out.done, samples_ch)

    emit:
    qc = qc_raw.out.fastqc_results
}
