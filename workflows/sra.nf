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

process prefetch_sra{
    conda "$FETCHENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)              "LOGS/$CONDITION/FETCH/Prefetch_SRA.log"
        else null
    }

    input:
    val collect
    val reads

    output:
    path "*.sra", emit: sra

    script:
        fn = reads+".sra"
        """
        prefetch $reads -o $fn &> prefetch.log
        """
}

process download_sra{
    conda "$FETCHENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)                "FASTQ/$CONDITION/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf(".log") >0)              "LOGS/$CONDITION/FETCH/SRA.log"
        else null
    }

    input:
    path sras

    output:
    path "*fastq.gz", emit: fq

    script:
    if (PAIRED == 'paired'){        
        """
        fasterq-dump -e $THREADS $FETCHPARAMS --split-files $sras &> sra.log ; rename 's/(.sra)*_([1|2])/_R$2/' *.fastq; for i in *.fastq;do pigz -p $THREADS \$i;done
        """
    }
    else{
        """
        fasterq-dump -e $THREADS $FETCHPARAMS $sras &> sra.log ; rename 's/(.sra)*_([1|2])/_R$2/' *.fastq ; for i in *.fastq;do pigz -p $THREADS \$i;done
        """
    }
}

workflow FETCH{
    take: collection

    main:
    //SAMPLE CHANNELS
    samples_ch = Channel.of(SHORTSAMPLES)

    collect_tofetch(collection.collect())
    prefetch_sra(collect_tofetch.out.done, samples_ch)
    download_sra(prefetch_sra.out.sra)

    emit:
    fetched = download_sra.out.fq
}
