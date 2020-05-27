FQSAMPLES = null
TRSAMPLES = null

if (PAIRED == 'paired'){
    FR1 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R1.fastq.gz"
    }
    FR2 = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_R2.fastq.gz"
    }
    FQSAMPLES = FR1+FR2
    FQSAMPLES.sort()

    TR1 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
    }
    TR2 = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
    }
    TRSAMPLES = TR1+TR2
    TRSAMPLES.sort()

}else{
    FQSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".fastq.gz"
    }
    FQSAMPLES.sort()
    TRSAMPLES=LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
    }
    TRSAMPLES.sort()
}

process qc_raw{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
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

process qc_trimmed{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
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
    path "*.{zip,html}", emit: trfastqc_results

    script:
    """
    fastqc --quiet -t $THREADS --noextract -f fastq $read
    """
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_raw{
    input:
    path results
    output:
    path "../QC/Multi/TRIMMED_RAW/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u > !{};done
    '''
}

//collecting list of processed file for multiqc, not implemented yet
process collect_qc_trimmed{
    input:
    path results
    output:
    path "../QC/Multi/TRIMMED_RAW/$CONDITION/qclist.txt", emit: collect_fastqc
    shell:
    '''
    for i in !{results};do echo $(dirname ${i}) >> tmp;done; cat tmp |sort -u > !{};done
    '''
}

process multiqc_trimmed_raw{
    conda "nextsnakes/envs/$TOOLENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1
    publishDir "${workflow.workDir}../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/TRIMMED_RAW/$CONDITION/$filename"
        else if (filename.indexOf("html") > 0)    "QC/Multi/TRIMMED_RAW/$CONDITION/$filename"
        else null
    }

    input:
    path qcs
    path trimmed
    output:
    path "*.{zip,html}", emit: multiqc_results

    script:
    """
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z ${workflow.workDir}/QC/FASTQC/$CONDITION/.
    """
}

workflow fastqc_trim{
    samples_ch = Channel.from(FQSAMPLES)
    trsamples_ch = Channel.from(TRSAMPLES)

    main:
    qc_raw(samples_ch)
    qc_trimmed(trsamples_ch)

    multiqc_raw(qc_raw.out.fastqc_results, qc_trimmed.out.trfastqc_results)
    emit:
    multiqc_trimmed_raw.out.multiqc_results
    //collect_qc_raw()
    //multiqc_raw(collect_qc_raw.out.collect_fastqc)
}
