COUNTENV = get_always('COUNTINGENV')
COUNTBIN = get_always('COUNTINGBIN')
COUNTREF = get_always('COUNTINGREF')
COUNTREFDIR = "${workflow.workDir}/../"+get_always('COUNTINGREFDIR')
COUNTANNO = get_always('COUNTINGANNO')
COUNTPREFIX = get_always('COUNTINGPREFIX') ?: COUNTBIN.split(' ')[0]

COUNTPARAMS = get_always('alevinfry_params_COUNT') ?: ''
PERMITPARAMS = get_always('alevinfry_params_PERMIT') ?: '--unfiltered-pl'

//COUNTING PROCESSES

process generate_t2g{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus 1
	cache 'lenient'

    input:
    path anno

    output:
    path "t2g.tsv", emit: t2g

    script:
    """
    zcat $anno | awk 'BEGIN{FS="\\t";OFS="\\t"} \$3=="transcript"{match(\$9,/transcript_id "[^"]+"/); t=substr(\$9,RSTART+15,RLENGTH-16); match(\$9,/gene_id "[^"]+"/); g=substr(\$9,RSTART+9,RLENGTH-10); print t,g}' > t2g.tsv
    """
}

process alevinfry_quant{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${COMBO}/${CONDITION}/COUNTING/alevinfry/${file(filename).getName()}"
        else                                    "COUNTS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path reads

    output:
    path "*_afquant", emit: counts
    path "*.log", emit: logs

    script:
    t2g = reads[0]
    raddir = reads[1]
    fn = file(raddir).getSimpleName().replaceAll(/\Q_map\E/,"")
    qd = fn+"_afquant"
    lf = "alevinfry_"+fn+".log"
    """
    $COUNTBIN generate-permit-list -i $raddir -d fw $PERMITPARAMS -o $qd &> $lf; $COUNTBIN collate -i $qd -r $raddir -t ${task.cpus} &>> $lf; $COUNTBIN quant -i $qd -m $t2g -t ${task.cpus} -r cr-like $COUNTPARAMS -o $qd --use-mtx &>> $lf
    """
}

workflow COUNTING{
    take: collection

    main:

    RADSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"/map.rad"
    }

    radsamples_ch = Channel.fromPath(RADSAMPLES.sort())

    annofile = Channel.fromPath(COUNTANNO)
    generate_t2g(annofile)
    alevinfry_quant(generate_t2g.out.t2g.combine(radsamples_ch.map{it -> file(it).getParent()}))

    emit:
    counts = alevinfry_quant.out.counts
    logs = alevinfry_quant.out.logs
}
