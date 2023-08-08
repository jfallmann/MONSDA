MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
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
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "hisat2.idx")            "$MAPIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/MAPPING/hisat2_index.log"
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
    zcat $gen > tmp.fa && mkdir -p $MAPUIDXNAME && $indexbin $IDXPARAMS -p ${task.cpus} tmp.fa $MAPUIDXNAME/$MAPPREFIX  &> index.log && ln -fs $MAPUIDXNAME hisat2.idx
    """

}

process hisat2_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)     "UNMAPPED/${COMBO}/${CONDITION}/"+"${file(filename).getName()}"
        else if (filename.indexOf(".summary") >0)            "MAPPED/${COMBO}/${CONDITION}/"+"${filename.replaceAll(/trimmed./,"")}"
        else if (filename.indexOf(".log") >0)               "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getName()}"
        else null
    }

    input:
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
    
    idx = reads[0]
    if (PAIRED == 'paired'){
        if (STRANDED == 'fr'){
            stranded = '--rna-strandness FR --fr'
        }else if (STRANDED == 'rf'){
            stranded = '--rna-strandness RF --rf'
        }else if (STRANDED == 'ff'){
            stranded = '--ff'
        }else{
            stranded = ''
        }

        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        pf = fn+"_mapped.sam"
        ufo = fn+"_R1_unmapped.fastq.gz"
        uft = fn+"_R2_unmapped.fastq.gz"
        lf = "hisat2_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS $stranded -p ${task.cpus} -x ${idx}/${MAPPREFIX} -1 $r1 -2 $r2 --un-conc-gz ${fn}.unmapped -S $pf &> $lf && gzip $pf && touch $ufo $uft && rename 's/.unmapped.([1|2]).gz/_R\$1_unmapped.fastq.gz/' ${fn}.unmapped.*.gz  &>> $lf
        """
    }else{
        if (STRANDED == 'fr'){
            stranded = '--rna-strandness F'
        }else if (STRANDED == 'rf'){
            stranded = '--rna-strandness R'
        }else{
            stranded = ''
        }
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "hisat2_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS $stranded -p ${task.cpus} -x ${idx}/${MAPPREFIX} -U $read --un-gz $uf -S $pf &> $lf && gzip $pf && touch $uf &>> $lf
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
        hisat2_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        hisat2_idx(genomefile)
        hisat2_mapping(hisat2_idx.out.htidx.combine(collection))
    }

    emit:
    mapped = hisat2_mapping.out.maps
    logs = hisat2_mapping.out.logs
}
