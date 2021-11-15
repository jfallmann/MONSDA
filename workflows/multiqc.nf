QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_MULTI') ?: ''

process collect_multi{
    input:
    path qclog
    path maplog
    path unique

    output:
    path "*mqccollect.txt", emit: done

    script:
    """
    echo "$qclog, $maplog, $unique Collection successful!" > mqccollect.txt
    """
}


process mqc{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else "QC/Multi/$COMBO$CONDITION/${file(filename).getName()}"
    }

    input:
    path others
    path samples

    output:
    path "*.zip", emit: mqc
    path "*.html", emit: html

    script:
    """
    touch $samples; export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o \${PWD} .
    """
}

process collect_dummy{
    input:
    path check

    output:
    path "mqc.txt", emit: done

    script:
    """
    echo "$check MQC successful!" > mqc.txt
    """
}

workflow MULTIQC{
    take:
    otherqcs
    maplogs
    unique

    main:

    //SAMPLE CHANNELS
    SAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"**_fastqc.zip"
    }
    SAMPLES.sort()
    samples_ch = Channel.fromPath(SAMPLES).ifEmpty([])
    
    MAPLOG = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../LOGS/$COMBO/MAPPING/**.log"
    }
    MAPLOG.sort()

    logs_ch = Channel.fromPath(MAPLOG).ifEmpty([])

    toqc_ch = samples_ch.mix(logs_ch)

    collect_multi(otherqcs.collect(), maplogs.collect(), unique.collect())
    mqc(collect_multi.out.done.collect(), toqc_ch)
    collect_dummy(mqc.out.mqc.collect())

    emit:
    mqcres = mqc.out.mqc
}
