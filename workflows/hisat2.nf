MAPENV=params.MAPPINGENV ?: null
MAPBIN=params.MAPPINGBIN ?: null

MAPIDX=params.MAPPINGIDX ?: null
MAPREF=params.MAPPINGREF ?: null
MAPGEN=params.MAPPINGGEN ?: null

IDXPARAMS = params.hisat2_params_0 ?: ''
MAPPARAMS = params.hisat2_params_1 ?: ''

//PROCESSES
process hisat_idx{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".ht2") > 0)        "$MAPIDX"+filename.replaceFirst(/tmp.idx/, '')
        else if (filename.indexOf(".idx") > 0)   "$MAPIDX"
        else null
    }

    input:
    val collect
    path reads
    path genome

    output:
    path "*.idx", emit: idx
    path "*.ht2", emit: htidx

    script:
    indexbin=MAPBIN.split(' ')[0]+'-build'
    """
    zcat $genome > tmp.fa && $indexbin $IDXPARAMS -p $THREADS tmp.fa tmp.idx && touch tmp.idx
    """

}

process hisat_mapping{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)          "UNMAPPED/$CONDITION/$filename"
        else if (filename.indexOf(".sam.gz") >0)        "MAPPED/$CONDITION/$filename"
        else if (filename.indexOf(".log") >0)           "MAPPED/$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path idxdummy
    path idx
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*fastq.gz", emit: unmapped
    path "*.log", emit: logs

    script:
    fn = file(reads[0]).getSimpleName()
    pf = fn+".mapped.sam"
    uf = fn+'.fastq.gz'
    if (STRANDED == 'fr'){
        stranded = '--rna-strandness F'
    }else if (STRANDED == 'rf'){
        stranded = '--rna-strandness R'
    }else{
        stranded = ''
    }

    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x $idx -1 $r1 -2 $r2 -S $pf --un-conc-gz $uf &> Log.log && gzip $pf && touch unmapped.fastq.gz
        """
    }else{
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x $idx -U $reads -S $pf --un-conc-gz $uf &> Log.log && gzip $pf && touch unmapped.fastq.gz
        """
    }
}

workflow MAPPING{
    take: samples_ch

    main:
    collect_results(samples_ch.collect())
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).merge(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    checkidx = file(MAPIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPIDX)
        hisat_mapping(collect_results.out.done, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        hisat_idx(collect_results.out.done, trimmed_samples_ch, genomefile)
        idxfile = Channel.fromPath(MAPIDX)
        hisat_mapping(collect_results.out.done, hisat_idx.out.idx, idxfile, trimmed_samples_ch)
    }


    emit:
    mapped  = hisat_mapping.out.maps
}
