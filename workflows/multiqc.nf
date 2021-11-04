QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_MULTI') ?: ''

process collect_multi{
    input:
    path check
    val checker

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}


process mqc{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else "QC/Multi/$COMBO$CONDITION/${file(filename).getName()}"
    }

    input:
    path others
    path samples

    output:
    path "*.{zip,html}", emit: mqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o \${PWD} .
    """
}

workflow MULTIQC{
    take:
    otherqcs
    maplogs

    main:

    //SAMPLE CHANNELS
    RSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_fastqc.zip"
    }
    RSAMPLES.sort()
    samples_ch = Channel.fromPath(RSAMPLES, followLinks: true).ifEmpty([])

    if (RUNDEDUP == 'enabled'){
        DSAMPLES=LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_dedup_fastqc.zip"
        }
        DSAMPLES.sort()
        dsamples_ch = Channel.fromPath(DSAMPLES, followLinks: true).ifEmpty([])
        samples_ch.join(dsamples_ch)
    }

    TSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_trimmed_fastqc.zip"
    }
    TSAMPLES.sort()
    tsamples_ch = Channel.fromPath(TSAMPLES, followLinks: true).ifEmpty([])

    MSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_mapped_sorted_fastqc.zip"
    }
    MSAMPLES.sort()

    USAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_mapped_sorted_unique_fastqc.zip"
    }
    USAMPLES.sort()

    MAPLOG = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../LOGS/$COMBO/MAPPING/*.log"
    }
    MAPLOG.sort()

    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true).ifEmpty([])
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true).ifEmpty([])
    logs_ch = Channel.fromPath(MAPLOG, followLinks: true).ifEmpty([])

    toqc_ch = samples_ch.concat(tsamples_ch, msamples_ch, usamples_ch, logs_ch)

    collect_multi(otherqcs.collect(), maplogs.collect())
    mqc(collect_multi.out.done.collect(), toqc_ch)

    emit:
    mqcres = mqc.out.mqc_results
}
