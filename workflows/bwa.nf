MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('bwa_params_INDEX') ?: ''
MAPPARAMS = get_always('bwa_params_MAP') ?: ''

IDXBIN = MAPBIN.split('_')[0]
MAPBIN = MAPBIN.replace('_', ' ')

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
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "bwa.idx")                          "$MAPIDX"
        else if (filename.indexOf("Log.out") > 0)           "LOGS/${COMBO}/${CONDITION}/bwa_index.log"
        else                                                "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "bwa.idx", emit: idx
    path "bwa_*", emit: bwidx

    script:
    gen =  genome.getName()
    """
    mkdir -p $MAPUIDXNAME && $IDXBIN index -p $MAPUIDXNAME/$MAPPREFIX $IDXPARAMS $gen &> Log.out && ln -fs $MAPUIDXNAME bwa.idx
    """

}

process bwa_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)   "UNMAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        //else if (filename.indexOf(".sam.gz") >0)          "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName().replaceAll(/_trimmed/,"")}"
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
    idx = reads[0]
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf1 = fn+"_R1_unmapped.fastq.gz"
        uf2 = fn+"_R2_unmapped.fastq.gz"
        lf = "bwa_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t ${task.cpus} ${idx}/${MAPPREFIX} $r1 $r2  2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> {log} &>/dev/null && touch $uf1 $uf2 2>> $lf &> /dev/null
        """
    }else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "bwa_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t ${task.cpus} ${idx}/${MAPPREFIX} $read  2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf
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
        bwa_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        bwa_idx(genomefile)
        bwa_mapping(bwa_idx.out.bwidx.combine(collection))
    }

    emit:
    mapped  = bwa_mapping.out.maps
    logs = bwa_mapping.out.logs
}
