MAPENV=params.MAPPINGENV ?: null
MAPBIN=params.MAPPINGBIN ?: null

MAPIDX=params.MAPPINGIDX ?: null
MAPUIDX=params.MAPPINGUIDX ?: null
MAPUIDXNAME=params.MAPPINGUIDXNAME ?: null
MAPREF=params.MAPPINGREF ?: null
MAPREFDIR=params.MAPPINGREFDIR ?: null
MAPANNO=params.MAPPINGANNO ?: null
MAPPREFIX=params.MAPPINGPREFIX ?: '.'

MAPUIDX.replace('.idx','')

IDXPARAMS = params.star_params_0 ?: ''
MAPPARAMS = params.star_params_1 ?: ''


//MAPPING PROCESSES

process collect_tomap{
    input:
    path check

    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}

process star_idx{
    conda "${workflow.workDir}/../NextSnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename =~ /SA/)                         "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename == "Genome")                "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".txt") > 0)        "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".tab") > 0)        "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf("Log.out") >0)      "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".idx") > 0)        "$MAPIDX"
        else null
    }

    input:
    val collect
    path reads
    path genome
    path anno

    output:
    path "*SA*", emit: idx
    path "*Log.out", emit: idxlog
    path "*.txt", emit: txts
    path "*.tab", emit: tabs
    path "*.idx", emit: tmpidx
    path "*Genome*", emit: idxgen

    script:
    gen =  genome.getName()
    an  = anno.getName()

    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && $MAPBIN $IDXPARAMS --runThreadN $THREADS --runMode genomeGenerate --outTmpDir STARTMP --genomeDir . --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && touch $MAPUIDXNAME && ln -s $MAPUIDXNAME star.idx
    """

}

process star_mapping{
    conda "${workflow.workDir}/../NextSnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("Unmapped.out") > 0)       "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/"\${PREFIX}"/,".").replaceAll(/Unmapped.out.*.gz/,"fastq.gz")}"
        else if (filename.indexOf(".sam.gz") >0)     "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/"\${PREFIX}"/,".").replaceAll(/trimmed.Aligned.out"/,"mapped")}"
        else if (filename.indexOf(".out") >0)        "MAPPED/$COMBO$CONDITION/$filename"
        else if (filename.indexOf(".tab") >0)        "MAPPED/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path idx
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*.out", emit: log
    path "*.tab", emit: sjtab
    path "*Unmapped.out*gz", includeInputs:false, emit: unmapped

    script:
    fn = file(reads[0]).getSimpleName()
    pf = fn+MAPPREFIX
    of = fn+MAPPREFIX+'Aligned.out.sam'

    """
    $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir ${workflow.workDir}/../$MAPUIDX --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix $pf --outReadsUnmapped Fastx && gzip $of && gzip *Unmapped.out*
    """
}

workflow MAPPING{
    take: collection

    main:
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
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    checkidx = file(MAPUIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPIDX)
        collect_tomap(collection.collect())
        star_mapping(collect_tomap.out.done, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        annofile = Channel.fromPath(MAPANNO)
        collect_tomap(collection.collect())
        star_idx(collect_tomap.out.done, trimmed_samples_ch, genomefile, annofile)
        star_mapping(collect_tomap.out.done, star_idx.out.idx, trimmed_samples_ch)
    }


    emit:
    mapped  = star_mapping.out.maps
}
