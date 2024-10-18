COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTIDX = get_always('COUNTINGIDX')
COUNTUIDX = get_always('COUNTINGUIDX')
COUNTUIDXNAME = get_always('COUNTINGUIDXNAME')+'.idx'
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTDECOY = get_always('COUNTINDECOY')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]

IDXPARAMS = get_always('kallisto_params_INDEX') ?: ''
COUNTPARAMS = get_always('kallisto_params_COUNT') ?: ''

//COUNTING PROCESSES

process trim{
    //conda "$TOOLENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_trimmed.fastq.gz") > 0)     "TRIMMED_FASTQ/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/${COMBO}/${CONDITION}/Trimming_report.txt"
        else null
    }

    input:
    path reads

    output:
    path "*trimmed.fastq.gz" , emit: trim
    path "Trimming_report.txt", emit: rep

    script:
    if (PAIRED == 'paired'){
        rs = reads[1..2].sort()
        r1 = rs[1]
        r2 = rs[2]
        a="Trimming_report.txt"
        b=file(r1).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        c=file(r2).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        """
        ln -sf $r1 $b ; ln -sf $r2 $c; echo "simulated $r1 $r2 trimming" > $a
        """
    }else{
        a="Trimming_report.txt"
        b=file(reads).getName().replace(".fastq.gz", "_trimmed.fastq.gz")
        """
        ln -sf $reads $b ; echo "simulated $reads trimming" > $a
        """
    }
}

process kallisto_idx{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "kallisto.idx")            "$COUNTIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/COUNTING/kallisto_index.log"
        else                                        "$COUNTUIDX"
    }

    input:
    path genome

    output:
    path "kallisto.idx", emit: idx

    script:    
    gen =  genome.getName()
    if (${COUNTINGDECOY}){
        decoy = "-d "+"${COUNTINGDECOY}" 
    }else{
        decoy = ''
    }
    """
    $COUNTBIN index $IDXPARAMS $decoy -t ${task.cpus} -i $COUNTUIDXNAME $gen &> index.log && ln -fs $COUNTUIDXNAME kallisto.idx
    """

}

process kallisto_quant{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${SCOMBO}/kallisto/${CONDITION}/COUNTING/${file(filename).getName()}"
        else                                    "COUNTS/${SCOMBO}/kallisto/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path reads

    output:
    path "*.gz", includeInputs:false, emit: counts
    path "*.log", emit: logs
    path "*.json", emit: json
    path "*.h5", emit: h5
    path "*", includeInputs:false, emit: rest

    script:

    idx = reads[0]
    if (PAIRED == 'paired'){
        if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '--fr-stranded'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
            stranded = '--rf-stranded'
        }else{
            stranded = ''
        }
        rs = reads[1..2].sort { a,b -> a[0] <=> b[0] == 0 ? (a[1..-1] as int) <=> (b[1..-1] as int) : a[0] <=> b[0] }
        r1 = rs[0]
        r2 = rs[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        lf = "kallisto_"+fn+".log"
        of = fn+"/abundance.tsv"
        oz = fn+"/abundance.tsv.gz"
        ol = fn+"_counts.gz"
        """
        $COUNTBIN quant -t ${task.cpus} -i $idx $stranded $COUNTPARAMS -o $fn $r1 $r2 &>> $lf && gzip $of && ln -fs $oz $ol && ln -fs ${fn}/run_info.json ${fn}_run_info.json && ln -fs ${fn}/abundance.h5 ${fn}_abundance.h5
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
        lf = "kallisto_"+fn+".log"
        of = fn+"/abundance.tsv"
        oz = fn+"/abundance.tsv.gz"
        ol = fn+"_counts.gz"
        """
         $COUNTBIN quant -t ${task.cpus} -i $idx $stranded $COUNTPARAMS -o $fn --single $read &>> $lf && gzip $of && ln -fs $oz $ol && ln -fs ${fn}/run_info.json ${fn}_run_info.json && ln -fs ${fn}/abundance.h5 ${fn}_abundance.h5
        """
    }
}

workflow COUNTING{
    take: collection

    main:
   
    checkidx = file(COUNTUIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(COUNTUIDX)
        if (PAIRED == 'paired'){
            kallisto_quant(idxfile.combine(samples_ch.collate(2)))
        } else{
            kallisto_quant(idxfile.combine(samples_ch.collate(1)))
        }        
    }
    else{
        genomefile = Channel.fromPath(COUNTREF)
        kallisto_idx(genomefile)
        if (PAIRED == 'paired'){
            kallisto_quant(kallisto_idx.out.idx.combine(samples_ch.collate(2)))
        } else{
            kallisto_quant(kallisto_idx.out.idx.combine(samples_ch.collate(1)))
        }
    }

    emit:
    counts = kallisto_quant.out.counts
    logs = kallisto_quant.out.logs
}
