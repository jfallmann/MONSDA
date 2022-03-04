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

IDXPARAMS = get_always('segemehl3_params_INDEX') ?: ''
MAPPARAMS = get_always('segemehl3_params_MAP') ?: ''


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
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "segemehl3.idx")                  "$MAPIDX"
        else if (filename.indexOf(".log") >0)             "LOGS/$COMBO$CONDITION/MAPPING/${file(filename).getName()}"
        else                                              "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "*.idx", emit: idx
    path "$MAPUIDXNAME", emit: uidx

    script:
    gen =  genome.getName()
    """
    $MAPBIN $IDXPARAMS --threads $THREADS -d $gen -x $MAPUIDXNAME &> index.log && ln -s $MAPUIDXNAME segemehl3.idx
    """

}

process segemehl3_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)   "UNMAPPED/$COMBO$CONDITION/"+"${file(filename).getSimpleName().replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
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
    genome = reads[0]
    idxfile = reads[1]
    gen =  genome.getName()
    idx = idxfile.getName()

    if (PAIRED == 'paired'){
        r1 = reads[2]
        r2 = reads[3]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $gen -q $r1 -p $r2 -o $pf -u $uf &> $lf && touch $uf && gzip *.sam
        """
    }else{
        fn = file(reads[2]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        read = reads[2]
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads $THREADS -i $idx -d $gen -q $read -o $pf -u $uf  &> $lf && touch $uf && gzip *.sam
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
        segemehl3_mapping(genomefile.combine(idxfile.combine(collection)))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        segemehl3_idx(genomefile)
        segemehl3_mapping(genomefile.combine(segemehl3_idx.out.idx.combine(collection)))
    }

    emit:
    mapped  = segemehl3_mapping.out.maps
    logs = segemehl3_mapping.out.logs
}
