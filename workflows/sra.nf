FETCHENV=get_always('FETCHENV')
FETCHBIN=get_always('FETCHBIN')

FETCHPARAMS = get_always('sra_params_DOWNLOAD') ?: ''


//FETCH PROCESSES

process collect_tofetch{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process fetchsra{
    conda "$FETCHENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)                "FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf(".log") >0)              "FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    val collect
    val reads

    output:
    path "*fastq.gz", emit: fq

    script:
    if (PAIRED == 'paired'){        
        """
        fasterq-dump -e {threads} $FETCHPARAMS --split-files $reads &> sra.log && rename 's/_1/_R1/' *.fastq && rename 's/_2/_R2/' *.fastq && pigz -p {threads} *.fastq
        """
    }
    else{
        """
        fasterq-dump -e {threads} $FETCHPARAMS $reads &> sra.log && rename 's/_1/_R1/' *.fastq && rename 's/_2/_R2/' *.fastq && pigz -p {threads} *.fastq
        """
    }
}

workflow FETCH{
    take: collection

    main:
    //SAMPLE CHANNELS
    samples_ch = Channel.of(SAMPLES)

    collect_tofetch(collection.collect())
    fetchsra(collect_tofetch.out.done, samples_ch)

    emit:
    fetched = fetchsra.out.fq
}
