MAPENV=params.MAPPINGENV ?: null
MAPBIN=params.MAPPINGBIN ?: null

MAPIDX=params.MAPPINGIDX ?: null
MAPREF=params.MAPPINGREF ?: null
MAPGEN=params.MAPPINGGEN ?: null
MAPANNO=params.MAPPINGANNO ?: null

IDXPARAMS = params.star_params_0 ?: ''
MAPPARAMS = params.star_params_1 ?: ''

//PROCESSES
process star_idx{
    conda "${workflow.workDir}/../nextsnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("SA") > 0)               "$MAPIDX"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".txt") > 0)        "$MAPIDX"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".tab") > 0)        "$MAPIDX"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf("Log.out") >0)      "$MAPIDX"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".idx") > 0)        "$MAPIDX"
        else null
    }

    input:
    val collect
    path reads
    path genome
    path anno

    output:
    path "SA*", emit: idx
    path "*Log.out", emit: idxlog
    path "*.txt", emit: txts
    path "*.tab", emit: tabs
    path "*.idx", emit: tmpidx

    script:
    gen =  genome.getName()
    an  = anno.getName()

    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && $MAPBIN $IDXPARAMS --runThreadN $THREADS --runMode genomeGenerate --outTmpDir STARTMP --genomeDir . --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && touch tmp.idx
    """

}

process star_mapping{
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
    val collect
    path idx
    path reads

    output:
    path "*Aligned.out.sam.gz", emit: maps
    path "*Log.out", emit: log
    path "*fastq.gz", includeInputs:false, emit: unmapped

    script:
    fn = file(reads[0]).getSimpleName()
    pf = fn+"."
    of = fn+'.Aligned.out.sam'

    """
    $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir $MAPGEN --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix $pf --outReadsUnmapped Fastx && gzip $of
    """
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
        star_mapping(collect_results.out.done, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        annofile = Channel.fromPath(MAPANNO)
        star_idx(collect_results.out.done, trimmed_samples_ch, genomefile, annofile)
        star_mapping(collect_results.out.done, star_idx.out.idx, trimmed_samples_ch)
    }


    emit:
    mapped  = star_mapping.out.maps
}
