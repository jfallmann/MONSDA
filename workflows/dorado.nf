CALLERENV = get_always('BASECALLENV')
CALLERBIN = get_always('BASECALLBIN')

CALLERPARAMS = get_always('dorado_params_CALLER') ?: ''
MODELPARAMS = get_always('dorado_params_MODEL') ?: ''

//CALLERS PROCESSES

process dorado{
    conda "$CALLERENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") > 0)        "LOGS/BASECALL/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path f5

    output:
    path "*.fastq.gz", emit: fastq
    path "*.log", emit: log

    script:
    fn = file(f5).getSimpleName()
    oc = fn+".fastq.gz"
    ol = fn+".log"
    sortmem = '30%'
    
    """
    $CALLERBIN download --model $MODELPARAMS &> $ol && $CALLERBIN basecaller $CALLERPARAMS $MODELPARAMS . 2>> $ol 1> tmp.bam && samtools view -h tmp.bam|samtools fastq -n - | pigz 1> $oc 2>> $ol && rm -rf tmp.bam
    """
}

workflow BASECALL{ 
    take: collection

    main:

    P5SAMPLES = SAMPLES.collect{
        element -> return "${workflow.workDir}/../RAW/"+element+"*.pod5"
    }

    p5samples_ch = Channel.fromPath(P5SAMPLES.sort())  
    
    dorado(p5samples_ch.collate(1))

    emit:
    fastq = dorado.out.fastq
    logs = dorado.out.log
}