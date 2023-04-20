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
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)             "LOGS/${COMBO}/${CONDITION}/star_index.log"
        else if (filename.indexOf(".idx") > 0)           "$MAPIDX"
        else if (filename == "$MAPUIDXNAME")             "$MAPUIDX"
        else                                             "$MAPUIDX/${filename}"
    }

    input:
    path genome
    path anno

    output:
    path "$MAPUIDXNAME", emit: idx
    path "*.out", emit: idxlog
    path "*.idx", emit: tmpidx

    script:
    gen =  genome.getName()
    an  = anno.getName()

    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p $MAPUIDXNAME && $MAPBIN $IDXPARAMS --runThreadN $THREADS --runMode genomeGenerate --outTmpDir STARTMP --genomeDir $MAPUIDXNAME --genomeFastaFiles tmp.fa --sjdbGTFfile tmp_anno && touch $MAPUIDXNAME && ln -s $MAPUIDXNAME star.idx && rm -f tmp.fa tmp_anno && ln -s $MAPUIDXNAME/* . && cat $MAPUIDXNAME/*Log.out > Log.out
    """
}

process star_mapping{
    conda "$MAPENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped") > 0)       "UNMAPPED/${COMBO}/${CONDITION}/"+"${file(filename).getName()}"
        //else if (filename.indexOf(".sam.gz") >0)     "MAPPED/${COMBO}/${CONDITION}/"+"${filename.replaceAll(/\Q.Aligned.out.sam.gz\E/,"")}_mapped.sam.gz"
        else if (filename.indexOf(".out") >0)        "LOGS/${COMBO}/${CONDITION}/MAPPING/star_"+"${filename.replaceAll(/\Q.out\E/,"")}.log"
        else if (filename.indexOf(".tab") >0)        "MAPPED/${COMBO}/${CONDITION}/"+"${filename}"
        else null
    }

    input:
    path reads

    output:
    path "*_mapped.sam.gz", emit: maps
    path "*.out", emit: logs
    path "*.tab", emit: sjtab
    path "*_unmapped.fastq.gz", includeInputs:false, emit: unmapped

    script:
    idx = reads[0]
    idxdir = idx.toRealPath()
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        a = "Trimming_report.txt"
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")+"."
        of = fn+'Aligned.out.sam'
        gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
        """
        $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2 --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed -r 's/.mate([1|2]).gz/_R\\1.gz/'| sed -r 's/\\.Unmapped.out_R([1|2]).gz/_R\\1_unmapped.fastq.gz/')"; done && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done
        """
    }
    else{
        if (PAIRED != 'singlecell'){
            read = reads[1]
            fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")+"."
            of = fn+'Aligned.out.sam'
            gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
            """
            $MAPBIN $MAPPARAMS --runThreadN $THREADS --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $read --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed -r 's/\\.Unmapped.out.mate1.gz/_unmapped.fastq.gz/')"; done && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done
            """
        }
        else{
            if (STRANDED == 'fr'){
                stranded = '--soloStrand Forward'
            }else if (STRANDED == 'rf'){
                stranded = '--soloStrand Reverse'
            }else{
                stranded = '--soloStrand Unstranded'
            }
            read = reads[1]
            fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
            umis = "${workflow.workDir}/../FASTQ/${CONDITION}/"+file(reads[1]).getSimpleName().replaceAll(/\QR2_trimmed\E/,"R1.fastq.gz")
            of = fn+'.Aligned.sortedByCoord.out.bam'
            gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")
            uf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_unmapped.fastq.gz")
            od = "${workflow.workDir}/../MAPPED/${COMBO}/${CONDITION}"

            """
            $MAPBIN --soloType CB_UMI_Simple $MAPPARAMS $stranded --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN $THREADS --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $read $umis --outFileNamePrefix ${fn}. --outReadsUnmapped Fastx && samtools view -h ${of} | gzip > $gf && rm -f $of ; paste <(cat ${fn}.Unmapped.out.mate1 | paste - - - -) <(cat ${fn}.Unmapped.out.mate2| paste - - - -) |tr \"\\t\" \"\\n\"| gzip > ${uf} && for f in *.Log.final.out; do mv "\$f" "\$(echo "\$f" | sed 's/.Log.final.out/.out/')"; done && mkdir -p $od && rsync -auv ${fn}.Solo.out $od
            """
        }
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
