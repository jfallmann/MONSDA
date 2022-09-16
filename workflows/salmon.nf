COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTIDX = get_always('COUNTINGIDX')
COUNTUIDX = get_always('COUNTINGUIDX')
COUNTUIDXNAME = get_always('COUNTINGUIDXNAME')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]
COUNTUIDX.replace('.idx','')

IDXPARAMS = get_always('salmon_params_INDEX') ?: ''
COUNTPARAMS = get_always('salmon_params_COUNT') ?: ''

//COUNTING PROCESSES

process salmon_idx{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "salmon.idx")            "$COUNTUIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/$COMBO$CONDITION/COUNTING/salmon_index.log"
    }

    input:
    path genome

    output:
    path "salmon.idx", emit: idx

    script:    
    gen =  genome.getName()
    """
    $COUNTBIN index $IDXPARAMS -p $THREADS -t $gen -i $COUNTUIDXNAME &> index.log && ln -fs $COUNTUIDXNAME salmon.idx
    """

}

process salmon_quant{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sf.gz") >0)            "COUNTS/$COMBO$CONDITION/"+"${filename.replaceAll(/trimmed./,"")}"
        else if (filename.indexOf(".log") >0)               "LOGS/$COMBO$CONDITION/COUNTING/${file(filename).getName()}"
        else null
    }

    input:
    path reads

    output:
    path "*.sf.gz", emit: counts
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
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        lf = "salmon_"+fn+".log"
        of = fn+"/COUNT.sf"
        oz = fn+"/COUNT.sf.gz"
        ol = fn+"_COUNT.sf.gz"
        """
        $COUNTBIN $COUNTPARAMS COUNT -p $THREADS -i $idx $stranded -o $fn -1 $r1 -2 $r2 &>> $lf && gzip $of && ln -sf $oz $ol
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
        of = fn+"/COUNT.sf"
        oz = fn+"/COUNT.sf.gz"
        ol = fn+"_COUNT.sf.gz"
        """
        $COUNTBIN $COUNTPARAMS COUNT -p $THREADS -i $idx $stranded -o $fn -r $read &>> $lf && gzip $of && ln -sf $oz $ol
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
        salmon_quant(idxfile.combine(samples_ch))
    }
    else{
        genomefile = Channel.fromPath(COUNTREF)
        salmon_idx(genomefile)
        salmon_quant(salmon_idx.out.idx.combine(samples_ch))
    }

    emit:
    counts = salmon_quant.out.counts
    logs = salmon_quant.out.logs
}
