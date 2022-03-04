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

IDXPARAMS = get_always('star_params_INDEX') ?: ''
MAPPARAMS = get_always('star_params_MAP') ?: ''

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

process star_idx{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)             "LOGS/$COMBO$CONDITION/star_index.log"
        else if (filename.indexOf(".idx") > 0)           "$MAPIDX"
        else if (filename == "$MAPUIDXNAME")             "$MAPUIDX"
        else                                             "$MAPUIDX/${filename}"
    }

    input:
    path genome
    path anno

    output:
    path "$MAPUIDXNAME", emit: idx
    path "*Log.out", emit: idxlog
    path "*.idx", emit: tmpidx

    script:
    gen =  genome.getName()
    an  = anno.getName()

    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p $MAPUIDXNAME && $MAPBIN $IDXPARAMS --runThreadN $THREADS --runMode genomeGenerate --outTmpDir STARTMP --genomeDir $MAPUIDXNAME --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && touch $MAPUIDXNAME && ln -s $MAPUIDXNAME star.idx && rm -f tmp.fa tmp_anno
    """

}

process star_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("Unmapped.out") > 0)       "UNMAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/\Q_trimmed.unmapped.out.gz\E/,"")}.fastq.gz"
        else if (filename.indexOf(".sam.gz") >0)     "MAPPED/$COMBO$CONDITION/"+"${filename.replaceAll(/\Q.Aligned.out.sam.gz\E/,"")}_mapped.sam.gz"
        else if (filename.indexOf(".out") >0)        "LOGS/$COMBO$CONDITION/MAPPING/star_"+"${filename.replaceAll(/\Q.out\E/,"")}.log"
        else if (filename.indexOf(".tab") >0)        "MAPPED/$COMBO$CONDITION/"+"${filename}"
        else null
    }

    input:
    path reads

    output:
    path "*.sam.gz", emit: maps
    path "*.out", emit: logs
    path "*.tab", emit: sjtab
    path "*Unmapped.out*gz", includeInputs:false, emit: unmapped

    script:
    idx = reads[0]
    idxdir = idx.toRealPath()
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        a = "Trimming_report.txt"
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")+"."
        of = fn+'Aligned.out.sam'
 
        """
        $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2 --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed 's/.mate[1|2].gz/.gz/')"; done && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done
        """
    }
    else{
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")+"."
        of = fn+'Aligned.out.sam'

        """
        $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $read --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed 's/.mate[1|2].gz/.gz/')"; done && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done
        """
    }
}

workflow MAPPING{
    take: collection

    main:
    checkidx = file(MAPUIDX)
    //collection.filter(~/.fastq.gz/)

    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        star_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        annofile = Channel.fromPath(MAPANNO)
        star_idx(genomefile, annofile)
        star_mapping(star_idx.out.idx.combine(collection))
    }


    emit:
    mapped  = star_mapping.out.maps
    logs = star_mapping.out.logs
}
