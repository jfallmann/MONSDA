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

IDXPARAMS = get_always('segemehl3_params_0') ?: ''
MAPPARAMS = get_always('segemehl3_params_1') ?: ''


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

process segemehl3_idx{
    conda "${workflow.workDir}/../NextSnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "$MAPUIDXNAME")                   "$MAPUIDX"
        else if (filename.indexOf(".idx") > 0)            "$MAPIDX"
        else null
    }

    input:
    val collect
    path reads
    path genome

    output:
    path "*.idx", emit: idx

    script:
    gen =  genome.getName()
    """
    $MAPBIN $IDXPARAMS --threads $THREADS -d $gen -x $MAPUIDXNAME && ln -s $MAPUIDXNAME segemehl3.idx
    """

}

process segemehl3_mapping{
    conda "${workflow.workDir}/../NextSnakes/envs/$MAPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
        saveAs: {filename ->
        if (filename.indexOf(".unmapped.fastq.gz") > 0)   "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_trimmed/,"")}"
        else if (filename.indexOf("Log.out") >0)          "MAPPED/$COMBO$CONDITION/$filename"
        else null
    }

    input:
    val collect
    path genome
    path idxfile
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*Log.out", emit: logs
    path "*fastq.gz", includeInputs:false, emit: unmapped

    script:
    fn = file(reads[0]).getSimpleName()
    pf = fn+".mapped.sam"
    uf = fn+".unmapped.fastq.gz"
    gen =  genome.getName()
    idx = idxfile.getName()

    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $gen -q $r1 -p $r2 -o $pf -u $uf &> Log.out && touch $uf && gzip *.sam
        """
    }else{
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $gen -q $reads -o $pf -u $uf  &> Log.out && touch $uf && gzip *.sam
        """
    }
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

    checkidx = file(MAPIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPIDX)
        genomefile = Channel.fromPath(MAPREF)
        collect_tomap(collection.collect())
        segemehl3_mapping(collect_tomap.out.done, genomefile, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        collect_tomap(collection.collect())
        segemehl3_idx(collect_tomap.out.done, trimmed_samples_ch, genomefile)
        segemehl3_mapping(collect_tomap.out.done, genomefile, segemehl3_idx.out.idx, trimmed_samples_ch)
    }


    emit:
    mapped  = segemehl3_mapping.out.maps
    logs = segemehl3_mapping.out.logs
}
