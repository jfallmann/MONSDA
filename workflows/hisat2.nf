MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX') ?: MAPBIN.split(' ')[0]
MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('hisat2_params_INDEX') ?: ''
MAPPARAMS = get_always('hisat2_params_MAP') ?: ''

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

process hisat2_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "hisat2.idx")            "$MAPIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/$COMBO$CONDITION/MAPPING/hisat2_index.log"
        else                                     "$MAPUIDX"
    }

    input:
    //val collect
    //path reads
    path genome

    output:
    path "hisat2.idx", emit: idx
    path "hisat2_*", emit: htidx

    script:
    indexbin=MAPBIN.split(' ')[0]+'-build'
    gen =  genome.getName()
    """
    zcat $gen > tmp.fa && mkdir -p $MAPUIDXNAME && $indexbin $IDXPARAMS -p $THREADS tmp.fa $MAPUIDXNAME/$MAPPREFIX  &> index.log && ln -fs $MAPUIDXNAME hisat2.idx
    """

}

process hisat2_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".unmapped.fastq.gz") > 0)     "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)            "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/trimmed./,"")}"
        else if (filename.indexOf(".log") >0)               "LOGS/$COMBO$CONDITION/MAPPING/hisat2.log"
        else null
    }

    input:
    //val collect
    path idx
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*fastq.gz", includeInputs:false, emit: unmapped
    path "*.log", emit: logs

    script:
    if (' ' in MAPBIN){
        mapbin = MAPBIN.split(' ')[1]
    } else {
        mapbin = MAPBIN
    }
    fn = file(reads[0]).getSimpleName()
    pf = fn+".mapped.sam"
    uf = fn+".unmapped.fastq.gz"

    if (STRANDED == 'fr'){
        stranded = '--rna-strandness F'
    }else if (STRANDED == 'rf'){
        stranded = '--rna-strandness R'
    }else{
        stranded = ''
    }

    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[0]
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x ${idx}/${MAPPREFIX} -1 $r1 -2 $r2 -S $pf --un-conc-gz $uf &> hisat_map.log && gzip *.sam && touch $uf
        """
    }else{
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x ${idx}/${MAPPREFIX} -U $reads -S $pf --un-conc-gz $uf &> hisat_map.log && gzip *.sam && touch $uf
        """
    }
}

workflow MAPPING{
    take: collection

    main:
   
    checkidx = file(MAPIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        if (PAIRED == 'paired'){
            hisat2_mapping(idxfile, collection.collate(2))
        }else{
            hisat2_mapping(idxfile, collection.collate(1))
        }
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        hisat2_idx(genomefile)
        if (PAIRED == 'paired'){
            hisat2_mapping(hisat2_idx.out.htidx, collection.collate(2))
        }else{
            hisat2_mapping(hisat2_idx.out.htidx, collection.collate(1))
        }
    }


    emit:
    mapped = hisat2_mapping.out.maps
    logs = hisat2_mapping.out.logs
}
