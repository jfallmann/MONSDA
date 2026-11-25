CALLERENV = get_always('BASECALLENV')
CALLERBIN = get_always('BASECALLBIN')

CALLERPARAMS = get_always('guppy_params_BASECALL') ?: ''
MODELPARAMS = get_always('guppy_params_MODEL') ?: ''

BASECALLMAP = new groovy.json.JsonSlurper().parse(new File(get_always('BASECALL_MANIFEST')))

//CALLERS PROCESSES

process guppy{
    conda "$CALLERENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$CALLERENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_summary.txt") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_telemetry.js") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") > 0)        "LOGS/${CONDITION}/BASECALL/guppy/${file(filename).getName()}"
    }

    input:
    tuple val(sample), path(f5, stageAs: 'inputs/signal????.fast5')

    output:
    path "*.fastq.gz", emit: fastq
    path "*_telemetry.js", emit: telemetry
    path "*_summary.txt", emit: summary
    path "*.log", emit: log

    script:
    fn = sample.split('/').last()
    oc = fn+".fastq.gz"
    ol = fn+".log"
    
    """
    mkdir -p TMP
    for f in inputs/*.fast5; do basename "\$f"; done > f5list
    $CALLERBIN $CALLERPARAMS -c $MODELPARAMS --compress_fastq -i inputs --input_file_list f5list -s TMP 2> $ol && cat TMP/pass/fastq_runid_*.fastq.gz > $oc && cat TMP/*.log >> $ol && mv -f TMP/sequencing_summary.txt ${fn}_summary.txt && mv -f TMP/sequencing_telemetry.js ${fn}_telemetry.js && rm -rf TMP f5list
    """
}

workflow BASECALL{ 
    take: collection

    main:

    SIGNALS = SAMPLES.collect{
        element -> return [element, BASECALLMAP[element].collect{ p -> return "${workflow.workDir}/../"+p }]
    }

    signals_ch = Channel.fromList(SIGNALS)
    
    guppy(signals_ch)

    emit:
    fastq = guppy.out.fastq
    logs = guppy.out.log
}