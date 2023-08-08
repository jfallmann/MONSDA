FETCHENV=get_always('FETCHENV')
FETCHBIN=get_always('FETCHBIN')

FETCHPARAMS = get_always('sra_params_PREFETCH') ?: ''
DOWNPARAMS = get_always('sra_params_DOWNLOAD') ?: ''


//FETCH PROCESSES

process prefetch_sra{
    conda "$FETCHENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)              "LOGS/$CONDITION/FETCH/Prefetch_SRA.log"
        else null
    }

    input:
    val reads

    output:
    path "*.sra", emit: sra

    script:
        fn = reads+".sra"
        """
        export NCBI_SETTINGS=\"$FETCHPARAMS\"
        prefetch $reads -o $fn &> prefetch.log
        """
}

process download_sra{
    conda "$FETCHENV"+".yaml"
    cpus THREADS
	cache 'lenient'
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
        export NCBI_SETTINGS=\"$FETCHPARAMS\"
        fasterq-dump -e ${task.cpus} $DOWNPARAMS --split-files $sras &> sra.log ; rename 's/(.sra)*_([1|2])/_R\$2/' *.fastq; for i in *.fastq;do pigz -p ${task.cpus} \$i;done
        """
    }
    else{
        """
        export NCBI_SETTINGS=\"$FETCHPARAMS\"
        fasterq-dump -e ${task.cpus} $DOWNPARAMS $sras &> sra.log ; rename 's/(.sra)*_([1|2])/_R\$2/' *.fastq ; for i in *.fastq;do pigz -p ${task.cpus} \$i;done
        """
    }
}

workflow FETCH{
    take: collection

    main:
    //SAMPLE CHANNELS
    samples_ch = Channel.of(SHORTSAMPLES)

    prefetch_sra(samples_ch)
    download_sra(prefetch_sra.out.sra)

    emit:
    fetched = download_sra.out.fq
}
