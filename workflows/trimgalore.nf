TRIMENV=get_always('TRIMMINGENV')
TRIMBIN=get_always('TRIMMINGBIN')

TRIMPARAMS = get_always('trimgalore_params_TRIM') ?: ''
cores = min(THREADS,4)
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
    cpus cores
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
        $TRIMBIN --cores $THREADS --paired --gzip $TRIMPARAMS $r1 $r2 &> trim.log && rename 's/_dedup//g' *.fq.gz && rename 's/_R([1|2])_val_([1|2]).fq.gz/_R\\1_trimmed.fastq.gz/g' *.fq.gz && rename 's/.fastq.gz_trimming/_trimming/g' *.txt
        """
    }
    else{
        """
        $TRIMBIN --cores $THREADS --gzip $TRIMPARAMS $reads &> trim.log && rename 's/_dedup//g' *.fq.gz && rename 's/.fq.gz/.fastq.gz/g' *.fq.gz && rename 's/.fastq.gz_trimming/_trimming/g' *.txt
        """
    }
}

workflow TRIMMING{
    take: collection

    main:
    if ( PREDEDUP == 'enabled' ){
        if (PAIRED == 'paired'){
            trim(collection.collate( 2 ))
        } else{
            trim(collection.collate( 1 ))
        }
    } else if ( collection.toList().contains('MONSDA.log') || collection.isEmpty()){
        //SAMPLE CHANNELS
        if (PAIRED == 'paired'){
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
            }                    
            collection = Channel.fromPath(SAMPLES).collate( 2 )
        }else{            
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
            }            
            collection = Channel.fromPath(SAMPLES).collate( 1 )
        }
        trim(collection)
    } else{
        if (PAIRED == 'paired'){
            trim(collection.collate( 2 ))
        } else{
            trim(collection.collate( 1 ))
        }
    }

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}
