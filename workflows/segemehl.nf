MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')+'.idx'
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')+'.idx'
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('segemehl_params_INDEX') ?: ''
MAPPARAMS = get_always('segemehl_params_MAP') ?: ''


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

process segemehl_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "segemehl.idx")                  "$MAPIDX"
        else if (filename.indexOf(".log") >0)             "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getName()}"
        else                                              "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "segemehl.idx", emit: idx
    path "$MAPUIDXNAME", emit: uidx

    script:
    gen =  genome.getName()
    """
    $MAPBIN $IDXPARAMS --threads ${task.cpus} -d $gen -x $MAPUIDXNAME &> index.log && ln -s $MAPUIDXNAME segemehl.idx
    """

}

process segemehl_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)   "UNMAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getName()}"
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
        pf = fn+"_mapped.sam.gz"
        uf1 = fn+"_R1_unmapped.fastq.gz"
        uf2 = fn+"_R2_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads ${task.cpus} -i $idx -d $gen -q $r1 -p $r2  2> $lf| tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> $lf &>/dev/null && touch $uf1 $uf2
        """
    }else{
        fn = file(reads[2]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        read = reads[2]
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads ${task.cpus} -i $idx -d $gen -q $read 2> $lf| tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf
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
        segemehl_mapping(genomefile.combine(idxfile.combine(collection)))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        segemehl_idx(genomefile)
        segemehl_mapping(genomefile.combine(segemehl_idx.out.idx.combine(collection)))
    }

    emit:
    mapped  = segemehl_mapping.out.maps
    logs = segemehl_mapping.out.logs
}
