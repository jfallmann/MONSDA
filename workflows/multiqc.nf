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


process multiqc{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    path others
    path samples
    path tsamples
    path msamples
    path usamples
    path logs

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -s 
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
    samples_ch = Channel.fromPath(RSAMPLES, followLinks: true)

    if (RUNDEDUP == 'enabled'){
        DSAMPLES=LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_dedup_fastqc.zip"
        }
        DSAMPLES.sort()
        dsamples_ch = Channel.fromPath(DSAMPLES, followLinks: true)
        samples_ch.join(dsamples_ch)
    }

    TSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_trimmed_fastqc.zip"
    }
    TSAMPLES.sort()
    tsamples_ch = Channel.fromPath(TSAMPLES, followLinks: true)

    MSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_mapped_sorted_fastqc.zip"
    }
    MSAMPLES.sort()

    USAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../QC/$COMBO"+element+"_mapped_sorted_unique_fastqc.zip"
    }
    USAMPLES.sort()

    MAPLOG = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/$COMBO"+element+".log"
    }
    MAPLOG.sort()

    msamples_ch = Channel.fromPath(MSAMPLES, followLinks: true)
    usamples_ch = Channel.fromPath(USAMPLES, followLinks: true)
    logs_ch = Channel.fromPath(MAPLOG, followLinks: true)

    collect_multi(otherqcs.collect(), maplogs.collect())
    multiqc(collect_multi.out.done.collect(), samples_ch, tsamples_ch, msamples_ch, usamples_ch, logs_ch)

    emit:
    mqcres = multiqc.out.multiqc_results
}
