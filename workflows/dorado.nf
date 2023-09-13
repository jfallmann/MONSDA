CALLERENV = get_always('BASECALLENV')
CALLERBIN = get_always('BASECALLBIN')

CALLERPARAMS = get_always('dorado_params_CALLER') ?: ''
MODELPARAMS = get_always('dorado_params_MODEL') ?: ''

//CALLERS PROCESSES

process guppy{
    conda "$CALLERENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_summary.txt") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_telemetry.js") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") > 0)        "LOGS/BASECALL/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path f5

    output:
    path ".fastq.gz", emit: fastq
    path "*_telemetry.js", emit: telemetry
    path "*_summary.txt", emit: summary
    path "*.log", emit: log

    script:
    fn = file(f5).getSimpleName()
    oc = fn+".fastq.gz"
    ol = fn+".log"
    sortmem = '30%'
    
    """
    mkdir -p TMP; ln -s *.pod5 TMP/. && $CALLERBIN download --model $MODELPARAMS &> $ol && $CALLERBIN basecaller $CALLERPARAMS TMP/ 2>> $ol 1> tmp.bam && samtools view -h tmp.bam|samtools fastq -n - | pigz > $oc && cat TMP/*.log >> $ol && mv -f TMP/sequencing_summary.txt . &&  mv -f TMP/sequencing_telemetry.js . && rm -rf TMP && rm -rf tmp.bam
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
    fastq = guppy.out.fastq
    logs = guppy.out.log
}