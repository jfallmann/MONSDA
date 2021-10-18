TRIMENV=get_always('TRIMMINGENV')
TRIMBIN=get_always('TRIMMINGBIN')

TRIMPARAMS = get_always('cutadapt_params_TRIM') ?: ''

//TRIMMING PROCESSES

process collect_totrim{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process trim{
    conda "$TRIMENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".fq.gz") > 0)                "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else if (filename.indexOf(".log") >0)              "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    val collect
    path reads

    output:
    path "*fq.gz", emit: trim
    path "*trimming_report.txt", emit: rep

    script:
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        """
        $TRIMBIN --cores $THREADS --paired $TRIMPARAMS $r1 $r2
        """
    }
    else{
        """
        $TRIMBIN --cores $THREADS $TRIMPARAMS $reads
        """
    }
}

workflow TRIMMING{
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
        RSAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
        }
        RSAMPLES.sort()
        samples_ch = Channel.fromPath(RSAMPLES)
    }

    collect_totrim(collection.collect())
    trim(collect_totrim.out.done, samples_ch)

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}