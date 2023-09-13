CALLERENV = get_always('BASECALLENV')
CALLERBIN = get_always('BASECALLBIN')

CALLERPARAMS = get_always('guppy_params_BASECALL') ?: ''
MODELPARAMS = get_always('guppy_params_MODEL') ?: ''

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
    mkdir -p TMP; echo \"${f5}\" > f5list && $CALLERBIN $CALLERPARAMS  -c $MODELPARAMS --compress_fastq -i . --input_file_list f5list -s TMP 2> $ol && cat TMP/pass/fastq_runid_*.fastq.gz > $oc && cat TMP/*.log >> $ol && mv -f TMP/sequencing_summary.txt . &&  mv -f TMP/sequencing_telemetry.js . && rm -rf TMP
    """
}

workflow BASECALL{ 
    take: collection

    main:

    F5SAMPLES = SAMPLES.collect{
        element -> return "${workflow.workDir}/../RAW/"+element+"*.fast5"
    }

    f5samples_ch = Channel.fromPath(F5SAMPLES.sort())  
    
    guppy(f5samples_ch.collate(1))

    emit:
    fastq = guppy.out.fastq
    logs = guppy.out.log
}