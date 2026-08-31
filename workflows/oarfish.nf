COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]

COUNTPARAMS = get_always('oarfish_params_COUNT') ?: ''

//COUNTING PROCESSES

process oarfish_quant{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${COMBO}/${CONDITION}/COUNTING/oarfish/${file(filename).getName()}"
        else                                    "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path bam

    output:
    path "*.gz", includeInputs:false, emit: counts
    path "*.log", emit: logs
    path "*", includeInputs:false, emit: rest

    script:
    fn = file(bam).getSimpleName().replaceAll(/\Q_mapped_sorted\E/,"")
    lf = "oarfish_"+fn+".log"
    oz = fn+"/quant.sf.gz"
    ol = fn+"_counts.gz"
    """
    mkdir -p $fn; samtools collate -@ ${task.cpus} -o namecollated.bam $bam &>> $lf; $COUNTBIN --threads ${task.cpus} --genome-alignments namecollated.bam --annotation $COUNTANNO --reference $COUNTREF $COUNTPARAMS --output $fn/oarfish &>> $lf && rm -f namecollated.bam && gzip -c $fn/oarfish.quant > $oz && ln -fs $oz $ol
    """
}

workflow COUNTING{
    take: collection

    main:

    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())

    oarfish_quant(mapsamples_ch.collate(1))

    emit:
    counts = oarfish_quant.out.counts
    logs = oarfish_quant.out.logs
}
