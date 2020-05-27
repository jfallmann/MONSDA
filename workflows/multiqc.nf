TOOLENV=params.QCENV ?: null
TOOLBIN=params.QCBIN ?: null

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

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/../QC/FASTQC/$CONDITION/.
    """
}

workflow MULTIQC{
    main:
    multiqc()
    emit:
    multiqc.out.multiqc_results
}
