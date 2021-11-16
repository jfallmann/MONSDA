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

IDXPARAMS = get_always('bwa_params_INDEX') ?: ''
MAPPARAMS = get_always('bwa_params_MAP') ?: ''


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

process bwa_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "bwa.idx")                  "$MAPIDX"
        else if (filename.indexOf("Log.out") >0)    "LOGS/$COMBO$CONDITION/bwa_index.log"
        else                                        "$MAPUIDX/${filename}"
    }

    input:
    //val collect
    //path reads
    path genome

    output:
    path "$MAPUIDXNAME", emit: uidx

    script:
    gen =  genome.getName()
    """
    touch $MAPUIDXNAME && $MAPBIN index -p $MAPPREFIX $IDXPARAMS $genome &> Log.out && ln -fs $MAPUIDXNAME bwa.idx
    """

}

process bwa_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf(".unmapped.fastq.gz") > 0)   "UNMAPPED/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/unmapped.fastq.gz/,"")}.fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_trimmed/,"")}"
        else if (filename.indexOf("Log.out") >0)          "LOGS/$COMBO$CONDITION/MAPPING/bwa.log"
        else null
    }

    input:
    //val collect
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
    uf = fn.replaceAll(/trimmed.fastq.gz/,"")+".unmapped.fastq.gz"
    gen =  genome.getName()
    idx = idxfile.getName()

    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        """
        $MAPBIN $MAPPARAMS --threads $THREADS $idx $r1 $r2|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> Log.out && touch $uf && gzip *.sam
        """
    }else{
        """
        $MAPBIN $MAPPARAMS --threads $THREADS $idx $reads|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> Log.out && touch $uf && gzip *.sam
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

    checkidx = file(MAPUIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        genomefile = Channel.fromPath(MAPREF)
        //collect_tomap(collection.collect())
        //bwa_mapping(collect_tomap.out.done, genomefile, idxfile, trimmed_samples_ch)
        bwa_mapping(genomefile, idxfile, collection.collect())
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        //collect_tomap(collection.collect())
        //bwa_idx(collect_tomap.out.done, trimmed_samples_ch, genomefile)
        //bwa_mapping(collect_tomap.out.done, genomefile, bwa_idx.out.idx, trimmed_samples_ch)
        bwa_idx(genomefile)
        bwa_mapping(genomefile, bwa_idx.out.uidx, collection.collect())
    }


    emit:
    mapped  = bwa_mapping.out.maps
    logs = bwa_mapping.out.logs
}
