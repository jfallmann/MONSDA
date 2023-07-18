//POST MAPPING PROCESSES

process sortsam{
    conda "samtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sam.gz") > 0)     "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.sam.gz"
        else null
    }

    input:
    path map

    output:
    path "*_mapped_sorted.sam.gz", includeInputs:false, emit: sam

    script:
    fn = file(map[0]).getSimpleName()
    sorted = fn+'_sorted.sam.gz'
    //No Maxthread in nextflow 
    sortmem = '30%'
    """
    set +o pipefail; mkdir -p TMP; samtools view -H $map|grep -P '^@HD' |pigz -p ${task.cpus} -f > tmphead; samtools view -H $map|grep -P '^@SQ'|sort -t\$'\\t' -k1,1 -k2,2V |pigz -p ${task.cpus} -f >> tmphead ; samtools view -H $map|grep -P '^@RG'|pigz -p ${task.cpus} -f >> tmphead ; samtools view -H $map|grep -P '^@PG'|pigz -p ${task.cpus} -f >> tmphead ; export LC_ALL=C;zcat $map | grep -v \"^@\"|sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k3,3V -k4,4n - |pigz -p ${task.cpus} -f > tmpfile; cat tmphead tmpfile > $sorted && rm -f tmphead tmpfile ${workflow.workDir}/../MAPPED/${COMBO}/${CONDITION}/$map
    """
}

process sam2bam{
    conda "samtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    path sam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: logs

    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+".bam"
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads ${task.cpus} && samtools index $bam 2> sam2bam.log
    """
}

process uniqsam{
    conda "samtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("unique.sam.gz") > 0)   "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.sam.gz"
        else if (filename.indexOf(".log") > 0)       "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    path sam

    output:
    path "*_unique.sam.gz", includeInputs:false, emit: sam
    path "*.log", emit: logs

    script:
    fn = file(sam[0]).getSimpleName()
    uniq = fn+'_unique.sam.gz'
    if (COMBO.indexOf('-bwa') > 0){
        bwa = 'bwa'
    } else{
        bwa = ''
    }
    """
    $BINS/Shells/UniqueSam_woPicard.sh $sam $uniq ${task.cpus} $bwa 2> uniq.log
    """
}

process sam2bamuniq{
    conda "samtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
        else null
    }

    input:
    path sam

    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: logs

    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+'.bam'
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads ${task.cpus} && samtools index $bam 2> uniquebam.log
    """
}


workflow POSTMAPPING{
    take: collection

    main:
  
    sortsam(collection)
    sam2bam(sortsam.out.sam)
    uniqsam(sortsam.out.sam)
    sam2bamuniq(uniqsam.out.sam)

    emit:
    postmap  = sam2bam.out.bam
    postbai  = sam2bam.out.bai
    postmapuni = sam2bamuniq.out.bam
    postunibai = sam2bamuniq.out.bai
}
