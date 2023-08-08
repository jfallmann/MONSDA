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
COUNTUIDX.replace('.idx','')

IDXPARAMS = get_always('salmon_params_INDEX') ?: ''
COUNTPARAMS = get_always('salmon_params_COUNT') ?: ''

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

process salmon_idx{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "salmon.idx")            "$COUNTUIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/COUNTING/salmon_index.log"
    }

    input:
    path genome

    output:
    path "salmon.idx", emit: idx

    script:    
    gen =  genome.getName()
    if (${COUNTINGDECOY}){
        decoy = "-d "+"${COUNTINGDECOY}" 
    }else{
        decoy = ''
    }
    """
    $COUNTBIN index $IDXPARAMS $decoy -p ${task.cpus} -t $gen -i $COUNTUIDXNAME &> index.log && ln -fs $COUNTUIDXNAME salmon.idx
    """

}

process salmon_quant{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${SCOMBO}/salmon/${CONDITION}/COUNTING/${file(filename).getName()}"
        else                                    "DTU/${SCOMBO}/salmon/${CONDITION}/"+"${filename.replaceAll(/trimmed./,"")}"
    }

    input:
    path reads

    output:
    path "*.gz", emit: counts
    path "*.log", emit: logs

    script:

    idx = reads[0]
    if (PAIRED == 'paired'){
        if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '-l ISF'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
            stranded = '-l ISR'
        }else{
            stranded = '-l IU'
        }
        rs = reads[1..2].sort { a,b -> a[0] <=> b[0] == 0 ? (a[1..-1] as int) <=> (b[1..-1] as int) : a[0] <=> b[0] }
        r1 = rs[0]
        r2 = rs[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        lf = "salmon_"+fn+".log"
        of = fn+"/quant.sf"
        oz = fn+"/quant.sf.gz"
        ol = fn+"_counts.gz"
        """
        $COUNTBIN $COUNTPARAMS quant -p ${task.cpus} -i $idx $stranded -o $fn -1 $r1 -2 $r2 &>> $lf && gzip $of && ln -fs $oz $ol
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
        lf = "salmon_"+fn+".log"
        of = fn+"/quant.sf"
        oz = fn+"/quant.sf.gz"
        ol = fn+"_counts.gz"
        """
        $COUNTBIN $COUNTPARAMS quant -p ${task.cpus} -i $idx $stranded -o $fn -r $read &>> $lf && gzip $of && ln -fs $oz $ol
        """
    }
}

workflow COUNTING{
    take: collection

    main:
   
    checkidx = file(COUNTIDX)
    collection.filter(~/.fastq.gz/)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(COUNTUIDX)
        if (PAIRED == 'paired'){
            salmon_quant(idxfile.combine(samples_ch.collate(2)))
        } else{
            salmon_quant(idxfile.combine(samples_ch.collate(1)))
        }        
    }
    else{
        genomefile = Channel.fromPath(COUNTREF)
        salmon_idx(genomefile)
        if (PAIRED == 'paired'){
            salmon_quant(salmon_idx.out.idx.combine(samples_ch.collate(2)))
        } else{
            salmon_quant(salmon_idx.out.idx.combine(samples_ch.collate(1)))
        }
    }

    emit:
    counts = salmon_quant.out.counts
    logs = salmon_quant.out.logs
}
