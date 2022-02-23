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

IDXPARAMS = get_always('minimap_params_INDEX') ?: ''
MAPPARAMS = get_always('minimap_params_MAP') ?: ''


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

process minimap_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "minimap.idx")                  "$MAPIDX"
        else if (filename.indexOf("Log.out") >0)          "LOGS/$COMBO$CONDITION/MAPPING/minimap_index.log"
        else                                        "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "*.idx", emit: idx
    path "$MAPUIDXNAME", emit: uidx

    script:
    gen =  genome.getName()
    """
    $MAPBIN -t $THREADS -d $MAPUIDXNAME $IDXPARAMS $genome &> index.log&& ln -fs $MAPUIDXNAME minimap.idx
    """

}

process minimap_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf(".unmapped.fastq.gz") > 0)   "UNMAPPED/$COMBO$CONDITION/${filename.replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_trimmed/,"")}"
        else if (filename.indexOf("Log.out") >0)          "LOGS/$COMBO$CONDITION/MAPPING/minimap.log"
        else null
    }

    input:
    path idxfile
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*Log.out", emit: logs
    path "*fastq.gz", includeInputs:false, emit: unmapped

    script:    
    idx = idxfile.getName()

    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+".mapped.sam"
        uf = fn+".unmapped.fastq.gz"
        """
        $MAPBIN $MAPPARAMS -t $THREADS $idx $r1 $r2|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> Log.out && touch $uf && gzip *.sam
        """
    }else{
        fn = file(reads).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+".mapped.sam"
        uf = fn+".unmapped.fastq.gz"
        """
        $MAPBIN $MAPPARAMS -t $THREADS $idx $reads|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> Log.out && touch $uf && gzip *.sam
        """
    }
}

workflow MAPPING{
    take: collection

    main:
 
    checkidx = file(MAPUIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        if (PAIRED == 'paired'){
            minimap_mapping(idxfile, collection)
        }else{
            minimap_mapping(idxfile, collection)
        }
        
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        minimap_idx(genomefile)
        if (PAIRED == 'paired'){
            minimap_mapping(minimap_idx.out.idx, collection)
        }else{
            minimap_mapping(minimap_idx.out.idx, collection)
        }
    }


    emit:
    mapped  = minimap_mapping.out.maps
    logs = minimap_mapping.out.logs
}
