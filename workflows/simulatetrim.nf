T1SAMPLES = null
T2SAMPLES = null

process trim{
    //conda "$TOOLENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("_trimmed.fastq.gz") > 0)     "TRIMMED_FASTQ/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/${COMBO}/${CONDITION}/Trimming_report.txt"
        else null
    }

    input:
    path reads

    output:
    path "*trimmed.fastq.gz" , emit: trim
    path "Trimming_report.txt", emit: rep

    script:
    if (PAIRED == 'paired' || PAIRED == 'singlecell'){
        r1 = reads[0]
        r2 = reads[1]
        a="Trimming_report.txt"
        b=file(r1).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        c=file(r2).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        """
        ln -sf $r1 $b ; ln -sf $r2 $c; echo "simulated $r1 $r2 trimming" > $a
        """
    }else{
        a="Trimming_report.txt"
        b=file(reads).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        """
        ln -sf $reads $b ; echo "simulated $reads trimming" > $a
        """
    }
}

workflow TRIMMING{
    take: 
    collection

    main:

    if ( PREDEDUP == 'enabled' ){
        trim(collection)
    } else if ( collection.toList().contains('MONSDA.log') || collection.isEmpty()){
        if (PAIRED == 'paired' || PAIRED == 'singlecell'){
            trim(samples_ch.collate(2))
        }
        else{
            trim(samples_ch.collate(1))
        }
    } else{
        trim(collection)
    }    

    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}