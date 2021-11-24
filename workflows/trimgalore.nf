TRIMENV=get_always('TRIMMINGENV')
TRIMBIN=get_always('TRIMMINGBIN')

TRIMPARAMS = get_always('trimgalore_params_TRIM') ?: ''

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
        if (filename.indexOf("_trimmed.fastq.gz") > 0)     "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed|_dedup/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/.fastq.gz/,"")}_trimming_report.txt"
        else if (filename.indexOf(".log") >0)              "LOGS/$COMBO$CONDITION/TRIMMING/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    //val collect
    path reads

    output:
    path "*_trimmed.fastq.gz", emit: trim
    path "*trimming_report.txt", emit: rep

    script:
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        """
        $TRIMBIN --cores $THREADS --paired --gzip $TRIMPARAMS $r1 $r2 &> trim.log && rename 's/_R([1|2])_val_([1|2]).fq.gz/_R\\1_trimmed.fastq.gz/g' *.fq.gz && rename 's/.fastq.gz_trimming/_trimming/g' *.txt
        """
    }
    else{
        """
        $TRIMBIN --cores $THREADS --gzip $TRIMPARAMS $reads &> trim.log && rename 's/.fq.gz/.fastq.gz/g' *.fq.gz && rename 's/.fastq.gz_trimming/_trimming/g' *.txt
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
                element -> return "${workflow.workDir}/../DEDUP_FASTQ/$COMBO"+element+"_dedup.*fastq.gz"
            }
        }
        else{
            SAMPLES = SAMPLES.collect{
            element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
            }
        }                 
    }

    if (collection.collect().contains('MONSDA.log') || collection.collect().isEmpty()){
        if (PAIRED == 'paired'){
            collection = Channel.fromPath(SAMPLES).collate( 2 )
        }
        else{
            collection = Channel.fromPath(SAMPLES).collate( 1 )
        }
    }

    trim(collection)

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}
