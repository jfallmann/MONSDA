TOOLENV=params.QCENV ?: null
TOOLBIN=params.QCBIN ?: null

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

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_raw{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_trimmed{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_map{
    input:
    path results
    output:
    path "QC/Multi/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u >> QC/Multi/!{$CONDITION}/qclist.txt;done
    '''
}

process multiqc{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1
    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path dummy

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z -s ${workflow.workDir}/../QC/FASTQC/${CONDITION}/.
    """
}

workflow MULTIQC{
    take: collection

    main:
    collect_multi(collection.collect())
    multiqc(collect_multi.out.done, collection.collect())

    emit:
    mqcres = multiqc.out.multiqc_results
}
