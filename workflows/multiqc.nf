TOOLENV=params.QCENV ?: null
TOOLBIN=params.QCBIN ?: null

MQCSAMPLES = null

if (PAIRED == 'paired'){
    R1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    R2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    MQCSAMPLES = R1+R2
    MQCSAMPLES.sort()

}else{
    MQCSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    MQCSAMPLES.sort()
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
    path dummy

    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/../QC/FASTQC/${CONDITION}/.
    """
}

workflow MULTIQC{
    take: dummy

    main:
    //mqcsamples_ch = Channel.from(MQCSAMPLES)
    multiqc(dummy)

    emit:
    mqcres = multiqc.out
}
