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

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "minimap.idx")                  "$MAPIDX"
        else if (filename.indexOf("Log.out") >0)          "LOGS/$COMBO$CONDITION/MAPPING/minimap_index.log"
        else                                        "$MAPUIDX"
    }

    input:
    //val collect
    //path reads
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
    uf = fn+".unmapped.fastq.gz"
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
 
    checkidx = file(MAPUIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        genomefile = Channel.fromPath(MAPREF)
        if (PAIRED == 'paired'){
            minimap_mapping(genomefile, idxfile, collection.collate(2))
        }else{
            minimap_mapping(genomefile, idxfile, collection.collate(1))
        }
        
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        minimap_idx(genomefile)
        if (PAIRED == 'paired'){
            minimap_mapping(genomefile, minimap_idx.out.idx, collection.collate(2))
        }else{
            minimap_mapping(genomefile, minimap_idx.out.idx, collection.collate(1))
        }
    }


    emit:
    mapped  = minimap_mapping.out.maps
    logs = minimap_mapping.out.logs
}
