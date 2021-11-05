MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('star_params_INDEX') ?: ''
MAPPARAMS = get_always('star_params_MAP') ?: ''

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
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename =~ /SA/)                         "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename == "Genome")                "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".txt") > 0)        "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf(".tab") > 0)        "$MAPUIDX"+"/"+"${filename.replaceAll(/star.idx/,"")}"
        else if (filename.indexOf("Log.out") >0)      "LOGS/$COMBO$CONDITION/star_index.log"
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
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf("Unmapped.out") > 0)       "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/\Q_trimmed.Unmapped.out.gz\E/,"")}.fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)     "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/\Qtrimmed.Aligned.out.sam.gz\E/,"")}mapped.sam.gz"
        else if (filename.indexOf(".out") >0)        "LOGS/$COMBO$CONDITION/MAPPING/star_"+"${filename.replaceAll(/\Q_trimmed\E/,"").replaceAll(/\Q.out\E/,"")}.log"
        else if (filename.indexOf(".tab") >0)        "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/\Q_trimmed\E/,"")}"
        else null
    }

    input:
    val collect
    path idx
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*.out", emit: logs
    path "*.tab", emit: sjtab
    path "*Unmapped.out*gz", includeInputs:false, emit: unmapped

    script:
    fn = file(reads[0]).getSimpleName()+'.'
    of = fn+'Aligned.out.sam'

    """
    $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir ${workflow.workDir}/../$MAPUIDX --readFilesCommand zcat --readFilesIn $reads --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed 's/.mate[1|2].gz/.gz/')"; done && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done
    """
}

workflow MAPPING{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    checkidx = file(MAPUIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
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
    logs = star_mapping.out.logs
}
