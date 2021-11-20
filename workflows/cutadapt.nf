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

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fq.gz") > 0)                "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed|_dedup/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}_trimming_report.txt"
        else if (filename.indexOf(".log") >0)              "LOGS/$COMBO$CONDITION/TRIMMING/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    //val collect
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
    //trim(samples_ch.collect())
    
    //collect_totrim(collection.collect())
    //trim(collect_totrim.out.done, samples_ch)
    collection.collect().join(samples_ch).unique().filter(!~/MONSDA.log/)
    trim(collection.collect())

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}
