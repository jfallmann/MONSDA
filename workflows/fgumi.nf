DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

EXTRACTPARAMS = get_always('fgumi_params_EXTRACT') ?: ''
DEDUPPARAMS = get_always('fgumi_params_DEDUP') ?: ''

process extract_fq{
    conda "$DEDUPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$DEDUPENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_dedup.fastq.gz") > 0)      "DEDUP_FASTQ/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf("log") > 0)             "LOGS/${COMBO}/${CONDITION}/DEDUP/dedup_extract.log"
        else null
    }

    input:
    path samples

    output:
    path "*_dedup.fastq.gz", emit: extract
    path "*_fgumi_extract.bam", emit: ubam
    path "ex.log", emit: logs

    script:
    if (PAIRED == 'paired'){
        r1 = samples[0]
        r2 = samples[1]
        sn = samples[0].getSimpleName().replace("_R1","")
        ubam = sn+"_fgumi_extract.bam"
        outf = samples[0].getSimpleName()+"_dedup.fastq.gz"
        outf2 = samples[1].getSimpleName()+"_dedup.fastq.gz"
        """
            mkdir -p tmp && $DEDUPBIN extract $EXTRACTPARAMS --inputs $r1 $r2 --sample $sn --library $sn --output $ubam > ex.log 2>&1 && samtools fastq -n -1 $outf -2 $outf2 -0 /dev/null -s /dev/null $ubam >> ex.log 2>&1
        """
    }
    else{
        r1 = samples[0]
        sn = samples[0].getSimpleName().replace(".fastq.gz","")
        ubam = sn+"_fgumi_extract.bam"
        outf = samples[0].getSimpleName()+"_dedup.fastq.gz"
        """
            mkdir -p tmp && $DEDUPBIN extract $EXTRACTPARAMS --inputs $r1 --sample $sn --library $sn --output $ubam > ex.log 2>&1 && samtools fastq -n $ubam | gzip -c > $outf && echo done >> ex.log
        """
    }
}

workflow DEDUPEXTRACT{
    take:
    collection

    main:
    //SAMPLE CHANNELS
    if ( PREDEDUP == 'enabled' ){
        if (PAIRED == 'paired'){
            extract_fq(samples_ch.collate( 2 ))
        } else{
            extract_fq(samples_ch.collate( 1 ))
        }
    }else{
        if (PAIRED == 'paired'){
            extract_fq(collection.collate( 2 ))
        } else{
            extract_fq(collection.collate( 1 ))
        }
    }

    emit:
    extract = extract_fq.out.extract
    ubam = extract_fq.out.ubam
    logs = extract_fq.out.logs
}

process dedup_bam{
    conda "$DEDUPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$DEDUPENV"
    cpus THREADS
    cache 'lenient'

    publishDir "${workflow.workDir}/../MAPPED/${COMBO}/${CONDITION}" , mode: 'link'

    input:
    path mapped_bam
    path ubam

    output:
    path "${mapped_bam.baseName}_dedup.bam", emit: dedup_bam
    path "${mapped_bam.baseName}_dedup.bam.bai", emit: dedup_idx
    path "dedup.log", emit: logs

    script:
    """
    mkdir -p tmp
    $DEDUPBIN zipper --unmapped $ubam --aligned $mapped_bam --output tmp/zippered.bam > dedup.log 2>&1
    $DEDUPBIN sort --order template-coordinate --input tmp/zippered.bam --output tmp/sorted.bam >> dedup.log 2>&1
    $DEDUPBIN dedup $DEDUPPARAMS --input tmp/sorted.bam --output ${mapped_bam.baseName}_dedup.bam >> dedup.log 2>&1
    samtools index ${mapped_bam.baseName}_dedup.bam >> dedup.log 2>&1
    rm $ubam
    """
}
