COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTIDX = get_always('COUNTINGIDX')
COUNTUIDX = get_always('COUNTINGUIDX')
COUNTUIDX = COUNTUIDX.replaceAll('.idx','')
COUNTUIDXNAME = get_always('COUNTINGUIDXNAME')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]

IDXPARAMS = get_always('simpleaf_params_INDEX') ?: '--ref-type spliced+unspliced'
COUNTPARAMS = get_always('simpleaf_params_COUNT') ?: '--chemistry 10xv3 --resolution cr-like --expected-ori fw --unfiltered-pl'

//COUNTING PROCESSES

process simpleaf_idx{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename == "simpleaf.idx")            "$COUNTIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/COUNTING/simpleaf/index.log"
        else                                      "$COUNTUIDX"
    }

    input:
    path genome
    path anno

    output:
    path "$COUNTUIDXNAME", emit: idx
    path "*.log", emit: idxlog

    script:
    gen =  genome.getName()
    an = anno.getName()
    """
    export ALEVIN_FRY_HOME=\$PWD/af_home; mkdir -p \$ALEVIN_FRY_HOME; $COUNTBIN set-paths &> index.log; zcat $gen > ref.fa && zcat $an > ref.gtf && $COUNTBIN index --fasta ref.fa --gtf ref.gtf -t ${task.cpus} $IDXPARAMS -o $COUNTUIDXNAME &>> index.log && rm -f ref.fa ref.gtf && ln -fs $COUNTUIDXNAME simpleaf.idx
    """

}

process simpleaf_quant{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${COMBO}/${CONDITION}/COUNTING/simpleaf/${file(filename).getName()}"
        else                                    "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path reads

    output:
    path "*_afquant", emit: counts
    path "*.log", emit: logs

    script:
    idx = reads[0]
    rs = reads[1..2].sort()
    r1 = rs[0]
    r2 = rs[1]
    fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
    qd = fn+"_afquant"
    lf = "simpleaf_"+fn+".log"
    """
    export ALEVIN_FRY_HOME=\$PWD/af_home; mkdir -p \$ALEVIN_FRY_HOME; $COUNTBIN set-paths &> $lf; $COUNTBIN quant --index $idx/index --reads1 $r1 --reads2 $r2 -t ${task.cpus} $COUNTPARAMS -o $qd &>> $lf
    """
}

workflow COUNTING{
    take: collection

    main:

    checkidx = file(COUNTUIDX)

    TRIMSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"_{R2,R1}_trimmed.fastq.gz"
    }

    trimsamples_ch = Channel.fromPath(TRIMSAMPLES.sort())

    if (checkidx.exists()){
        idxfile = Channel.fromPath(COUNTUIDX)
        simpleaf_quant(idxfile.combine(trimsamples_ch.collate(2)))
    }
    else{
        genomefile = Channel.fromPath(COUNTREF)
        annofile = Channel.fromPath(COUNTANNO)
        simpleaf_idx(genomefile, annofile)
        simpleaf_quant(simpleaf_idx.out.idx.combine(trimsamples_ch.collate(2)))
    }

    emit:
    counts = simpleaf_quant.out.counts
    logs = simpleaf_quant.out.logs
}
