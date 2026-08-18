
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
    container "oras://jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }

    input:
    path samples

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    OUT=\${PWD}
    LOGDIR="${workflow.workDir}/../LOGS"
    SCAN="\${LOGDIR}/${CONDITION}"
    COLLECT="\${SCAN}/MULTIQC/collect"
    VERSIONS="\${LOGDIR}/versions.txt"
    mkdir -p "\$COLLECT"

    for i in $samples; do
        cp -f "\$i" "\$COLLECT"/
    done

    MODS=""
    if [[ -f "\$VERSIONS" ]]; then
        MODS=\$(grep -v '^#' "\$VERSIONS" | cut -f3 | grep -vx '-' | sort -u | sed 's/^/-m /' | tr '\\n' ' ')
        cp -f "\$VERSIONS" "\$OUT"/
    fi
    export LC_ALL=C.UTF-8
    multiqc -f \$MODS -k json -z -s -o "\$OUT" "\$SCAN"
    """
}

workflow PREMULTIQC{
    take:
    otherqcs

    main:

    //SAMPLE CHANNELS
    premultiqc(otherqcs.collect())

    emit:
    mqcres = premultiqc.out.multiqc_results
}
