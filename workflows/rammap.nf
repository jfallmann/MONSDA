MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX = MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('rammap_params_INDEX') ?: ''
MAPPARAMS = get_always('rammap_params_MAP') ?: ''

//MAPPING PROCESSES

process rammap_idx{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "rammap.idx")                  "$MAPIDX"
        else if (filename.indexOf("index.log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/rammap_index.log"
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
    gzip -cdfq $gen > tmp_ref.fa && $MAPBIN -t ${task.cpus} -d $MAPUIDXNAME $IDXPARAMS tmp_ref.fa &> index.log && rm -f tmp_ref.fa && ln -fs $MAPUIDXNAME rammap.idx
    """

}

process rammap_mapping{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

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
    idxfile = reads[0]
    idx = idxfile.getName()
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf1 = fn+"_R1_unmapped.fastq.gz"
        uf2 = fn+"_R2_unmapped.fastq.gz"
        lf = "rammap_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t ${task.cpus} $idx $r1 $r2 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> $lf &>/dev/null && touch $uf1 $uf2 2>> $lf &> /dev/null
        """
    }else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "rammap_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -t ${task.cpus} $idx $read 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf
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
        rammap_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        rammap_idx(genomefile)        
        rammap_mapping(rammap_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = rammap_mapping.out.maps
    logs = rammap_mapping.out.logs
}
