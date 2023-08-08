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

IDXBIN = 'bwameth.py index-mem2' //MAPBIN.split('_')[0]
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

process bwameth_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)             "LOGS/${COMBO}/${CONDITION}/bwameth_index.log"
        else if (filename.indexOf(".idx") > 0)           "$MAPIDX"
        else                                             "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "$MAPUIDXNAME", emit: idx
    path "Log.out", emit: idxlog
    path "*.idx", emit: tmpidx

    script:
    gen =  genome.getName()
    genfa = MAPPREFIX+genome.getName().replace('.gz', '')
    """
    mkdir -p $MAPUIDXNAME && zcat $gen > $MAPUIDXNAME/$genfa && $IDXBIN $MAPUIDXNAME/$genfa $IDXPARAMS &> Log.out && ln -fs $MAPUIDXNAME/$genfa bwameth.idx
    """

}

process bwameth_mapping{
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
        lf = "bwameth_"+fn+".log"
        """
        touch ${idx}/${MAPPREFIX}*c2t*; $MAPBIN $MAPPARAMS --threads ${task.cpus} --reference ${idx}/${MAPPREFIX}*.fa $r1 $r2 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> $lf &>/dev/null && touch $uf1 $uf2 2>> $lf &> /dev/null
        """
    }else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "bwameth_"+fn+".log"
        """
        touch ${idx}/${MAPPREFIX}*c2t*; $MAPBIN $MAPPARAMS --threads ${task.cpus} --reference ${idx}/${MAPPREFIX}*.fa $read 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf
        """
    }
}

workflow MAPPING{
    take: collection

    main:
    
    checkidx = file(MAPUIDX)
    collection.filter(~/.fastq.gz/)

    if (checkidx.exists()){
        //idxfile = Channel.fromPath("${MAPUIDX}"+File.separator+"${MAPPREFIX}"+"*")
        idxfile = Channel.fromPath(MAPUIDX)
        bwameth_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        bwameth_idx(genomefile)
        bwameth_mapping(bwameth_idx.out.idx.combine(collection))
    }

    emit:
    mapped  = bwameth_mapping.out.maps
    logs = bwameth_mapping.out.logs
}
