QCENV=get_always('QCENV')
QCBIN=get_always('QCBIN')
QCPARAMS = get_always('fastqc_params_MULTI') ?: ''

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
    //path samples

    output:
    path "*.zip", emit: mqc
    path "*.html", emit: html

    script:
    """
    touch $others; export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -o \${PWD} .
    """
}

workflow MULTIQC{
    take:
    otherqcs
    
    main:
    
    mqc(otherqcs.collect())

    emit:
    mqcres = mqc.out.mqc
}
