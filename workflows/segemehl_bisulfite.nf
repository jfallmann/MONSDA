MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN').split('_')[0]
MAPIDX = get_always('MAPPINGIDX')
BISIDX = get_always('MAPPINGIDX2')
MAPUIDX = get_always('MAPPINGUIDX')+'.idx'
MAPUIDX2 = get_always('MAPPINGUIDX2')+'.idx2'
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')+'.idx'
MAPUIDX2NAME = get_always('MAPPINGUIDX2NAME')+'.idx2'
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('segemehlbisulfite_params_INDEX') ?: ''
MAPPARAMS = get_always('segemehlbisulfite_params_MAP') ?: ''
MAPENV = "${MAPENV}".replace('bisulfite', '')

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
        if (filename == "segemehlbisulfite.idx")                  "$MAPIDX"
        else if (filename == "segemehlbisulfite_bs.idx2")          "$BISIDX"
        else if (filename.indexOf(".log") >0)             "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getName()}"
        else if (filename == "$MAPUIDX2NAME")             "$MAPUIDX2"
        else                                              "$MAPUIDX"
        
    }

    input:
    path genome

    output:
    path "segemehlbisulfite.idx", emit: idx
    path "segemehlbisulfite_bs.idx2", emit: idx2
    path "$MAPUIDXNAME", emit: uidx
    path "$MAPUIDX2NAME", emit: uidx2

    script:
    gen =  genome.getName()
    """
    $MAPBIN $IDXPARAMS --threads ${task.cpus} -d $gen -x $MAPUIDXNAME -y $MAPUIDX2NAME &> index.log && ln -s $MAPUIDXNAME segemehlbisulfite.idx && ln -s $MAPUIDX2NAME segemehlbisulfite_bs.idx2
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
    idxfile2 = reads[2]
    gen =  genome.getName()
    idx = idxfile.getName()
    idx2 = idxfile2.getName()

    if (PAIRED == 'paired'){
        r1 = reads[3]
        r2 = reads[4]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf1 = fn+"_R1_unmapped.fastq.gz"
        uf2 = fn+"_R2_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads ${task.cpus} -i $idx -j $idx2 -d $gen -q $r1 -p $r2 -o tmp.sam 2> $lf && cat tmp.sam| tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> $lf &>/dev/null && touch $uf1 $uf2 && rm -f tmp.sam
        """
    }else{
        fn = file(reads[3]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        read = reads[3]
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "segemehl_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS --threads ${task.cpus} -i $idx -j $idx2 -d $gen -q $read -o tmp.sam 2> $lf && cat tmp.sam| tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf && rm -f tmp.sam
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
        idxfile2 = Channel.fromPath(MAPUIDX2)
        genomefile = Channel.fromPath(MAPREF)
        segemehl_mapping(genomefile.combine(idxfile.combine(idxfile2.combine(collection))))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        segemehl_idx(genomefile)
        segemehl_mapping(genomefile.combine(segemehl_idx.out.idx.combine(segemehl_idx.out.idx2.combine(collection))))
    }

    emit:
    mapped  = segemehl_mapping.out.maps
    logs = segemehl_mapping.out.logs
}
