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
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1
    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/$CONDITION/$filename"
        else null
    }

    input:
    path qcs
    path trimmed
    path mapped
    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/QC/FASTQC/$CONDITION/.
    """
}

workflow multiqc{
    main:
    multiqc_raw(qc_raw.out.fastqc_results, qc_trimmed.out.trfastqc_results, qc_mapped.out.mapfastqc_results)
    emit:
    multiqc.out.multiqc_results
    //collect_qc_raw()
    //multiqc_raw(collect_qc_raw.out.collect_fastqc)
}
