MQCENV=get_always('POSTQCENV')
MQCBIN=get_always('POSTQCBIN')
MQCPARAMS = get_always('rustqc_params_MULTI') ?: ''

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
    SCAN="\${LOGDIR}/${COMBO}/${CONDITION}"
    COLLECT="\${SCAN}/MULTIQC/collect"
    VERSIONS="\${LOGDIR}/versions.txt"
    EXTRA=""
    BASE_QC_DIR="${workflow.workDir}/../QC"
    COMBO_VAL="${COMBO}"
    CONDITION_VAL="${CONDITION}"
    mkdir -p "\$COLLECT"

    for i in $others; do
        cp -f "\$i" "\$COLLECT"/
    done

    # If the corresponding fastqc combo exists, include its output in the MultiQC report.
    FQ_COMBO="\${COMBO_VAL/rustqc/fastqc}"
    FQ_DIR="\${BASE_QC_DIR}/\${FQ_COMBO}/\${CONDITION_VAL}"
    if [[ "\$FQ_COMBO" != "\$COMBO_VAL" && -d "\$FQ_DIR" ]]; then
        EXTRA="\$FQ_DIR"
    fi

    MODS=""
    if [[ -f "\$VERSIONS" ]]; then
        MODS=\$(grep -v '^#' "\$VERSIONS" | cut -f3 | tr ',' '\\n' | grep -vx '-' | sort -u | sed 's/^/-m /' | tr '\\n' ' ')
        cp -f "\$VERSIONS" "\$OUT"/
    fi
    export LC_ALL=C.UTF-8
    multiqc -f $MQCPARAMS \$MODS -k json -z -s -o "\$OUT" "\$SCAN" \$EXTRA
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
