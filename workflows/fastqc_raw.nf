FQSAMPLES = null

if (PAIRED == 'paired'){
    R1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    R2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    FQSAMPLES = R1+R2
    FQSAMPLES.sort()

}else{
    FQSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    FQSAMPLES.sort()
}

process qc_raw{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$CONDITION/$filename"
        else null
    }

    input:
    path read

    output:
    path "*.{zip,html}", emit: fastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

process collect_qc_raw{
    input:
    path results
    output:
    path "../QC/Multi/RAW/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u > !{};done
    '''
}

process multiqc_raw{
    conda "${workflow.workDir}/../nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1
    publishDir "${workflow.workDir}../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/RAW/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/RAW/$CONDITION/$filename"
        else null
    }

    input:
    path qcs
    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/QC/FASTQC/$CONDITION/.
    """
}

workflow {
    samples_ch = Channel.from(FQSAMPLES)
    main:
    qc_raw(samples_ch)
    multiqc_raw(qc_raw.out.fastqc_results)
    emit:
    multiqc_raw.out.multiqc_results
    //collect_qc_raw()
    //multiqc_raw(collect_qc_raw.out.collect_fastqc)
}
