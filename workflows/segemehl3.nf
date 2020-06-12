MAPENV=params.MAPPINGENV ?: null
MAPBIN=params.MAPPINGBIN ?: null

MAPIDX=params.MAPPINGIDX ?: null
MAPREF=params.MAPPINGREF ?: null
MAPGEN=params.MAPPINGGEN ?: null

IDXPARAMS = params.segemehl3_params_0 ?: ''
MAPPARAMS = params.segemehl3_params_1 ?: ''

//PROCESSES
process sege_idx{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".idx") > 0)        "$MAPIDX"
        else null
    }

    input:
    //val collect
    path genome

    output:
    path "*.idx", emit: idx

    script:
    """
    $MAPBIN $IDXPARAMS --threads $THREADS -d $genome -x tmp.idx
    """

}

process sege_mapping{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)          "UNMAPPED/$CONDITION/$filename"
        else if (filename.indexOf(".sam.gz") >0)        "MAPPED/$CONDITION/$filename"
        else if (filename.indexOf("Log.out") >0)        "MAPPED/$CONDITION/$filename"
        else null
    }

    input:
    //val collect
    path genome
    path idx
    path reads

    output:
    path "TMP/STAROUT/*Aligned.out.sam.gz", emit: maps
    path "TMP/STAROUT/*Log.out", emit: maplog
    path "TMP/STAROUT/*fastq.gz", emit: unmapped

    script:
    fn = file(reads[0]).getSimpleName()
    pf = fn+".mapped.sam"
    uf = fn+'.fastq.gz'

    if (PAIRED == 'paired'){
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $genome -q $reads[0] -p $reads[1] -o $pf -u $uf
        """
    }else{
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $genome -q $reads -o $pf -u $uf
        """
    }
}

workflow MAPPING{
    take: samples_ch

    main:
    //collect_results(samples_ch.collect())
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

    checkidx = File(MAPIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPIDX)
        genomefile = Channel.fromPath(MAPREF)
        sege_mapping(genomefile, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        sege_idx(genomefile, samples_ch)
        sege_mapping(genomefile, sege_idx.out.idx, trimmed_samples_ch)
    }


    emit:
    mapped  = sege_mapping.out.maps
}
