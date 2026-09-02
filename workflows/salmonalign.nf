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
MAPYAML = MAPBIN

IDXPARAMS = get_always('salmonalign_params_INDEX') ?: ''
MAPPARAMS = get_always('salmonalign_params_MAP') ?: ''


//MAPPING PROCESSES

process salmonalign_idx{
    conda "$MAPYAML"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPYAML"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "salmonalign.idx")              "$MAPIDX"
        else if (filename.indexOf("index.log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/salmon/salmonalign_index.log"
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
    $MAPBIN index $IDXPARAMS -p ${task.cpus} -t $gen -i $MAPUIDXNAME &> index.log && ln -fs $MAPUIDXNAME salmonalign.idx
    """

}

process salmonalign_mapping{
    conda "$MAPYAML"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPYAML"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'

    publishDir "${workflow.workDir}/../" , mode: 'link',
        saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)   "UNMAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") >0)          "LOGS/${COMBO}/${CONDITION}/MAPPING/salmon/${file(filename).getName()}"
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
        if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '-l ISF'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
            stranded = '-l ISR'
        }else{
            stranded = '-l IU'
        }
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf1 = fn+"_R1_unmapped.fastq.gz"
        uf2 = fn+"_R2_unmapped.fastq.gz"
        lf = "salmonalign_"+fn+".log"
        """
        $MAPBIN quant -p ${task.cpus} -i $idx $stranded $MAPPARAMS --writeMappings=/dev/stdout -o $fn -1 $r1 -2 $r2 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools collate -u -O -|samtools fastq -n -c 6 -1 $uf1 -2 $uf2 ) 2>> $lf &>/dev/null && touch $uf1 $uf2 2>> $lf &> /dev/null
        """
    }else{
        if (STRANDED == 'fr' || STRANDED == 'SF'){
            stranded = '-l SF'
        }else if (STRANDED == 'rf' || STRANDED == 'SR'){
            stranded = '-l SR'
        }else{
            stranded = '-l U'
        }
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam.gz"
        uf = fn+"_unmapped.fastq.gz"
        lf = "salmonalign_"+fn+".log"
        """
        $MAPBIN quant -p ${task.cpus} -i $idx $stranded $MAPPARAMS --writeMappings=/dev/stdout -o $fn -r $read 2> $lf|tee >(samtools view -h -F 4 |gzip > $pf) >(samtools view -h -f 4 |samtools fastq -n - | pigz > $uf) 2>> $lf &> /dev/null && touch $uf
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
        salmonalign_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        salmonalign_idx(genomefile)        
        salmonalign_mapping(salmonalign_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = salmonalign_mapping.out.maps
    logs = salmonalign_mapping.out.logs
}
