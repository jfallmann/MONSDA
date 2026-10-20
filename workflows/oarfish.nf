BINS = get_always('BINS')
COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]

COUNTPARAMS = get_always('oarfish_params_COUNT') ?: ''

//COUNTING PROCESSES

process decompress_annotation_oarfish{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus 1
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("log") > 0)        "LOGS/${COMBO}/COUNTING/oarfish/decompress_annotation.log"
    }

    input:
    path anno

    output:
    path "annotation.gtf", emit: anno_plain
    path "log", emit: log

    script:
    """
    zcat -f $anno > annotation.gtf 2> log; if [ ! -s annotation.gtf ]; then echo "ERROR: decompressed annotation from $anno is empty, check source file and filesystem" >> log; exit 1; fi
    """
}

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
    path fls

    output:
    path "*.gz", includeInputs:false, emit: counts
    path "*.log", emit: logs
    path "*", includeInputs:false, emit: rest

    script:
    anno = fls[0]
    bam = fls[1]
    fn = file(bam).getSimpleName().replaceAll(/\Q_mapped_sorted\E/,"")
    lf = "oarfish_"+fn+".log"
    oz = fn+"/quant.sf.gz"
    ol = fn+"_counts.gz"
    """
    mkdir -p $fn; : > $lf; if [ ! -s $anno ]; then echo "ERROR: decompressed annotation $anno is empty, check source file and filesystem" >> $lf; exit 1; fi; samtools collate -@ ${task.cpus} -o namecollated.bam $bam &>> $lf; $COUNTBIN --threads ${task.cpus} --genome-alignments namecollated.bam --annotation $anno --genome-fasta $COUNTREF $COUNTPARAMS --output $fn/oarfish &>> $lf && rm -f namecollated.bam && python3 $BINS/Analysis/oarfish_to_salmon.py -i $fn/oarfish.quant -o $oz &>> $lf && ln -fs $oz $ol
    """
}

workflow COUNTING{
    take: collection

    main:

    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())
    annofile = Channel.fromPath(COUNTANNO)
    decompress_annotation_oarfish(annofile)
    annogtf = decompress_annotation_oarfish.out.anno_plain

    oarfish_quant(annogtf.combine(mapsamples_ch.collate(1)))

    emit:
    counts = oarfish_quant.out.counts
    logs = oarfish_quant.out.logs
}
