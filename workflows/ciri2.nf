CIRCENV = get_always('CIRCSENV')
CIRCBIN = get_always('CIRCSBIN')
CIRCREF = get_always('CIRCSREF')
CIRCREFDIR = "${workflow.workDir}/../"+get_always('CIRCSREFDIR')
CIRCANNO = get_always('CIRCSANNO')

CIRCPARAMS = get_always('ciri2_params_CIRC') ?: ''

//CIRCS PROCESSES

process ciri2{
    conda "$CIRCENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_circs") > 0)      "CIRCS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}"        
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}"
    }

    input:
    path fls

    output:
    path "*_circs", emit: circs
    path "log", emit: log

    script:
    ref = fls[0]
    anno = fls[1]
    reads = fls[2]        
    fn = file(reads).getSimpleName()
    oc = fn+"_circs"
    ol = fn+".log"
    sortmem = '30%'
    
    """
    set +o pipefail; export LC_ALL=C; if [[ -n \"\$(zcat ${reads} | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then mkdir -p TMP && zcat ${reads}|samtools sort -n -@ ${task.cpus} -u -O sam -T TMP > ${fn}_tmp.sam && zcat ${anno} > ${fn}_tmp.gtf && zcat ${ref} > ${fn}_tmp.fa && perl $CIRCBIN -I ${fn}_tmp.sam -O ${fn}_circs -F ${fn}_tmp.fa -T ${task.cpus} -A ${fn}_tmp.gtf -G log $CIRCPARAMS &>> log; else gzip < /dev/null > ${fn}_circs; echo \"File ${reads} empty\" >> log; fi; touch CIRIerror.log && cat CIRIerror.log >> {log} && echo '' > CIRIerror.log && touch ${fn}_circs
    """
}

workflow CIRCS{ 
    take: collection

    main:

    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"*_mapped_sorted.sam.gz"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())  
    annofile = Channel.fromPath(CIRCANNO)
    genomefile = Channel.fromPath(CIRCREF)

    ciri2(genomefile.combine(annofile.combine(mapsamples_ch.collate(1))))

    emit:
    circs = ciri2.out.circs
    logs = ciri2.out.log
}