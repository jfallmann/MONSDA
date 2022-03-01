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
        else if (filename.indexOf("index.log") >0)          "LOGS/$COMBO$CONDITION/MAPPING/minimap_index.log"
        else                                            "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "*.idx", emit: idx
    path "$MAPUIDXNAME", emit: uidx

    script:
    gen =  genome.getName()
    """
    $MAPBIN -t $THREADS -d $MAPUIDXNAME $IDXPARAMS $gen &> index.log && ln -fs $MAPUIDXNAME minimap.idx
    """

}

process minimap_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)   "UNMAPPED/$COMBO$CONDITION/${filename.replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)          "MAPPED/$COMBO$CONDITION/${file(filename).getSimpleName().replaceAll(/_trimmed/,"")}"
        else if (filename.indexOf(".log") >0)          "LOGS/$COMBO$CONDITION/MAPPING/${file(filename).getName()}"
        else null
    }

    input:
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*fastq.gz", includeInputs:false, emit: unmapped
    path "*.log", emit: logs

    script:    
    idxfile = reads[0]
    idx = idxfile.getName()
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "minimap_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t $THREADS $idx $r1 $r2|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> $lf && touch $uf && gzip *.sam
        """
    }else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "minimap_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t $THREADS $idx $reads|tee >(samtools view -h -F 4 > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 1>/dev/null 2&> $lf && touch $uf && gzip *.sam
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
        minimap_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        minimap_idx(genomefile)        
        minimap_mapping(minimap_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = minimap_mapping.out.maps
    logs = minimap_mapping.out.logs
}
