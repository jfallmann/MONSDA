MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX') ?: 'piscem_idx'

IDXPARAMS = get_always('piscem_params_INDEX') ?: ''
MAPPARAMS = get_always('piscem_params_MAP') ?: ''


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

process piscem_idx{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "piscem.idx")                  "$MAPIDX"
        else if (filename.indexOf("index.log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/piscem_index.log"
        else                                            "$MAPUIDXNAME"
    }

    input:
    path genome

    output:
    path "$MAPUIDXNAME", emit: idx
    path "*index.log", emit: idxlog

    script:
    gen =  genome.getName()
    """
    mkdir -p $MAPUIDXNAME && zcat $gen > ref.fa && $MAPBIN build -t ${task.cpus} $IDXPARAMS -s ref.fa -o $MAPUIDXNAME/$MAPPREFIX &> index.log && rm -f ref.fa && ln -fs $MAPUIDXNAME piscem.idx
    """

}

process piscem_mapping{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf(".log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getName()}"
        else                                       "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path reads

    output:
    path "*_map", emit: maps
    path "*.log", emit: logs

    script:
    idxdir = reads[0]
    r1 = reads[1]
    r2 = reads[2]
    fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
    of = fn+"_map"
    lf = "piscem_"+fn+".log"
    """
    $MAPBIN map-sc -t ${task.cpus} $MAPPARAMS -i $idxdir/$MAPPREFIX -1 $r1 -2 $r2 -o $of &> $lf
    """
}

workflow MAPPING{
    take: collection

    main:

    checkidx = file(MAPUIDX)
    collection.filter(~/.fastq.gz/)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        piscem_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        piscem_idx(genomefile)
        piscem_mapping(piscem_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = piscem_mapping.out.maps
    logs = piscem_mapping.out.logs
}
