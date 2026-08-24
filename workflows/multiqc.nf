MQCENV=get_always('PREQCENV')
MQCBIN=get_always('PREQCBIN')
MQCPARAMS = get_always('fastqc_params_MULTI') ?: ''

process mqc{
    conda "$MQCENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MQCENV"
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
    path "versions.txt", emit: versions, optional: true

    script:
    """
    touch $others
    OUT=\${PWD}
    LOGDIR="${workflow.workDir}/../LOGS"
    BASE_QC_DIR="${workflow.workDir}/../QC"
    COMBO_VAL="${COMBO}"
    CONDITION_VAL="${CONDITION}"
    VERSIONS="\${LOGDIR}/versions.txt"
    SCAN="\${LOGDIR}/\${COMBO_VAL}/\${CONDITION_VAL}"

    # All QC results of this combo are reported, no need to collect them first.
    QC_DIR="\${BASE_QC_DIR}/\${COMBO_VAL}/\${CONDITION_VAL}"
    if [[ -d "\$QC_DIR" ]]; then
        SCAN="\$SCAN \$QC_DIR"
    fi

    # If the corresponding fastqc combo exists, include its output in the MultiQC report.
    FQ_DIR="\${BASE_QC_DIR}/\${COMBO_VAL/rustqc/fastqc}/\${CONDITION_VAL}"
    if [[ "\$FQ_DIR" != "\$QC_DIR" && -d "\$FQ_DIR" ]]; then
        SCAN="\$SCAN \$FQ_DIR"
    fi

    MODS=""
    if [[ -f "\$VERSIONS" ]]; then
        MODS=\$(grep -v '^#' "\$VERSIONS" | cut -f3 | tr ',' '\\n' | grep -vx '-' | sort -u | sed 's/^/-m /' | tr '\\n' ' ')
        cp -f "\$VERSIONS" "\$OUT"/
    fi
    export LC_ALL=C.UTF-8
    multiqc -f $MQCPARAMS \$MODS -k json -z -s -o "\$OUT" \$SCAN
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
