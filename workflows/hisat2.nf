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
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".ht2") > 0)        "$MAPUIDX"+"/"+"${filename.replaceFirst(/tmp\.idx/, '')}"
        else if (filename.indexOf(".idx") > 0)   "$MAPIDX"
        else null
    }

    input:
    val collect
    path reads
    path genome

    output:
    path "*.idx", emit: idx
    path "*.ht2", emit: htidx

    script:
    indexbin=MAPBIN.split(' ')[0]+'-build'
    gen =  genome.getName()
    """
    zcat $gen > tmp.fa && $indexbin $IDXPARAMS -p $THREADS tmp.fa $MAPUIDXNAME/$indexbin && ln -fs $MAPUIDXNAME hisat2.idx
    """

}

process hisat2_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".unmapped.fastq.gz") > 0)     "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/unmapped.fastq.gz/,"")}fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)            "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/trimmed./,"")}"
        else if (filename.indexOf(".log") >0)               "MAPPED/$COMBO$CONDITION/"+"${filename}.log"
        else null
    }

    input:
    val collect
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
        r1 = reads[0]
        r2 = reads[1]
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x ${idx}/${MAPUIDXNAME} -1 $r1 -2 $r2 -S $pf --un-conc-gz $uf &> hisat_map.log && gzip *.sam && touch $uf
        """
    }else{
        """
        $MAPBIN $MAPPARAMS $stranded -p $THREADS -x ${idx}/${MAPUIDXNAME} -U $reads -S $pf --un-conc-gz $uf &> hisat_map.log && gzip *.sam && touch $uf
        """
    }
}

workflow MAPPING{
    take: collection

    main:
    //SAMPLE CHANNELS
    if (PAIRED == 'paired'){
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R1_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        T2SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_R2_trimmed.fastq.gz"
        }
        T2SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES).join(Channel.fromPath(T2SAMPLES))

    }else{
        T1SAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/$COMBO"+element+"_trimmed.fastq.gz"
        }
        T1SAMPLES.sort()
        trimmed_samples_ch = Channel.fromPath(T1SAMPLES)
    }

    checkidx = file(MAPIDX)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPIDX)
        collect_tomap(collection.collect())
        hisat2_mapping(collect_tomap.out.done, idxfile, trimmed_samples_ch)
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        collect_tomap(collection.collect())
        hisat2_idx(collect_tomap.out.done, trimmed_samples_ch, genomefile)
        hisat2_mapping(collect_tomap.out.done, hisat2_idx.out.htidx, trimmed_samples_ch)
    }


    emit:
    mapped = hisat2_mapping.out.maps
    logs = hisat2_mapping.out.logs
}
