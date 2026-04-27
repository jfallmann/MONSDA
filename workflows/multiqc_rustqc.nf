QCENV=get_always('POSTQCENV')
QCBIN=get_always('POSTQCBIN')
QCPARAMS = get_always('rustqc_params_MULTI') ?: ''

process mqc{
    conda "$QCENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path others, stageAs: 'mqc_input??/*'
    //path samples

    output:
    path "*.zip", emit: mqc
    path "*.html", emit: html

    script:
    """
    touch $others
    OUT=\${PWD}
    LIST=multiqc_inputs.txt
    TMP_LIST=multiqc_inputs_unique.txt
    BASE_QC_DIR="${workflow.workDir}/../QC"
    COMBO_VAL="${COMBO}"
    CONDITION_VAL="${CONDITION}"

    for i in $others; do
        dirname "\$i" >> "\$LIST"
    done

    # If the corresponding fastqc combo exists, include its output in the MultiQC report.
    FQ_COMBO="\${COMBO_VAL/rustqc/fastqc}"
    FQ_DIR="\${BASE_QC_DIR}/\${FQ_COMBO}/\${CONDITION_VAL}"
    if [[ -d "\$FQ_DIR" ]]; then
        echo "\$FQ_DIR" >> "\$LIST"
    fi

    sort -u "\$LIST" > "\$TMP_LIST"
    export LC_ALL=en_US.utf8
    export LC_ALL=C.UTF-8
    multiqc -f --exclude picard --exclude gatk -k json -z -s -o "\$OUT" -l "\$TMP_LIST"
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
