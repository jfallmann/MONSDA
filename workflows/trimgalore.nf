TRIMENV=get_always('TRIMMINGENV')
TRIMBIN=get_always('TRIMMINGBIN')

TRIMPARAMS = get_always('trimgalore_params_TRIM') ?: ''
//int cores = min(THREADS,4)
//TRIMMING PROCESSES

process trim{
    conda "$TRIMENV"+".yaml"
    cpus 4//cores
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
    take: 
    samples_ch    

    main:
    check = samples_ch.toList()
    if ( check.contains('MONSDA.log') || samples_ch.isEmpty()){
        //SAMPLE CHANNELS
        if (PAIRED == 'paired'){
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
            }                    
            collection = Channel.fromPath(SAMPLES)
            trim(collection.collate( 2 ))
        }else{            
            SAMPLES = SAMPLES.collect{
                element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
            }            
            collection = Channel.fromPath(SAMPLES)
            trim(collection.collate( 1 ))
        }
    } else if ( PREDEDUP == 'enabled' ){
        if (PAIRED == 'paired'){
            trim(samples_ch.collate(2))
        } else{
            trim(samples_ch.collate(1))
        }
    } else{
        collection = samples_ch
        if (PAIRED == 'paired'){
            trim(collection.collate(2))
        } else{
            trim(collection.collate(1))
        }
    }

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}
