MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX') ?: MAPBIN.split(' ')[0]
MAPUIDX = MAPUIDX.replace('.idx','')

IDXPARAMS = get_always('bowtie2_params_INDEX') ?: ''
MAPPARAMS = get_always('bowtie2_params_MAP') ?: ''

// NOTE: unlike hisat2, bowtie2 has no --rna-strandness/stranded-protocol
// flag, so MAPPINGSTRANDED is intentionally unused here.

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

process bowtie2_idx{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename == "bowtie2.idx")           "$MAPIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/MAPPING/bowtie2/index.log"
        else                                     "$MAPUIDX"
    }

    input:
    path genome

    output:
    path "bowtie2.idx", emit: idx
    path "bowtie2_*", emit: btidx

    script:
    indexbin=MAPBIN.split(' ')[0]+'-build'
    gen =  genome.getName()
    """
    zcat $gen > tmp.fa && mkdir -p $MAPUIDXNAME && $indexbin $IDXPARAMS --threads ${task.cpus} tmp.fa $MAPUIDXNAME/$MAPPREFIX  &> index.log && ln -fs $MAPUIDXNAME bowtie2.idx
    """

}

process bowtie2_mapping{
    conda "$MAPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped.fastq.gz") > 0)     "UNMAPPED/${COMBO}/${CONDITION}/"+"${file(filename).getName()}"
        else if (filename.indexOf(".log") >0)               "LOGS/${COMBO}/${CONDITION}/MAPPING/bowtie2/${file(filename).getName()}"
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
    if (MAPPINGPAIRED == 'paired'){
        rds = sort_reads(reads[1..2])
        r1 = rds[0]
        r2 = rds[1]
        fn = file(r1).getSimpleName().replaceAll(/_R1(_dedup)?_trimmed$/,"")
        pf = fn+"_mapped.sam"
        ufo = fn+"_R1_unmapped.fastq.gz"
        uft = fn+"_R2_unmapped.fastq.gz"
        lf = "bowtie2_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -p ${task.cpus} -x ${idx}/${MAPPREFIX} -1 $r1 -2 $r2 --un-conc-gz ${fn}.unmapped -S $pf &> $lf && gzip $pf && touch $ufo $uft && rename 's/.unmapped.([1|2]).gz/_R\$1_unmapped.fastq.gz/' ${fn}.unmapped.*.gz  &>> $lf
        """
    }else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/(_dedup)?_trimmed$/,"")
        pf = fn+"_mapped.sam"
        uf = fn+"_unmapped.fastq.gz"
        lf = "bowtie2_"+fn+".log"
        """
        $MAPBIN $MAPPARAMS -p ${task.cpus} -x ${idx}/${MAPPREFIX} -U $read --un-gz $uf -S $pf &> $lf && gzip $pf && touch $uf &>> $lf
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
        bowtie2_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        bowtie2_idx(genomefile)
        bowtie2_mapping(bowtie2_idx.out.btidx.combine(collection))
    }

    emit:
    mapped = bowtie2_mapping.out.maps
    logs = bowtie2_mapping.out.logs
}
