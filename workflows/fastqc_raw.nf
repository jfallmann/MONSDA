if (PAIRED == 'paired'){
    //log('Running paired end raw QC')

    Channel
        .from(SAMPLES+'_{R1,R2}.fastq.gz') //{ file -> file.name.replaceAll(/.bam|.bai$/,'') }
        //.set { samples_ch }
        .view()
}else{
    //log('Running single end raw QC')

    Channel
        .from(SAMPLES+'.fastq.gz')
        .view()
        //.set { samples_ch }
}

process qc_raw{
    conda 'nextsnakes/envs/qc.yaml'
    validExitStatus 0,1
    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/FASTQC/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$filename"
        else null
    }

    input:
    set prefix, file(reads) from samples_ch

    output:
    file "*.{zip,html}" into fastqc_results

    script:
    """
    echo ${prefix}
    fastqc $reads
    """
}

process multiqc_raw{
    conda 'nextsnakes/envs/qc.yaml'
    validExitStatus 0,1
    publishDir "${workflow.workDir}" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/RAW/$filename"
        else if (filename.indexOf("html") > 0)    "QC/FASTQC/$filename"
        else null
    }

    input:
    set sampleId, file(reads) from samples_ch
    output:

    output:
    file "*.{zip,html}" into fastqc_results

    script:
    """
    echo ${prefix}
    export LC_ALL=en_US.utf8; export LC_ALL=C.UTF-8; multiqc -f --exclude picard --exclude gatk -k json -z
    """
}
