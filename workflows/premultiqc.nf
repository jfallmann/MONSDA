
QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_MULTI') ?: ''

process collect_multi{
    input:
    path check
    
    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}


process premultiqc{
    conda "$QCENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$COMBO$CONDITION/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    path samples

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -s 
    """
}

workflow PREMULTIQC{
    take:
    otherqcs

    main:

    //SAMPLE CHANNELS
    multiqc(otherqcs.collect())

    emit:
    mqcres = premultiqc.out.multiqc_results
}
