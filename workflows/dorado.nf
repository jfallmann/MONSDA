CALLERENV = get_always('BASECALLENV')
CALLERBIN = get_always('BASECALLBIN')

CALLERPARAMS = get_always('dorado_params_BASECALL') ?: get_always('dorado_params_CALLER') ?: ''
MODELPARAMS = get_always('dorado_params_MODEL') ?: ''

BASECALLMAP = new groovy.json.JsonSlurper().parse(new File(get_always('BASECALL_MANIFEST')))

//CALLERS PROCESSES

process dorado{
    conda "$CALLERENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$CALLERENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".fastq.gz") > 0)      "FASTQ/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".log") > 0)        "LOGS/${CONDITION}/BASECALL/dorado/${file(filename).getName()}"
    }

    input:
    tuple val(sample), path(p5, stageAs: 'inputs/signal????.pod5')

    output:
    path "*.fastq.gz", emit: fastq
    path "*.log", emit: log

    script:
    fn = sample.split('/').last()
    oc = fn+".fastq.gz"
    ol = fn+".log"
    
    """
    set -o pipefail
    $CALLERBIN download --model $MODELPARAMS &> $ol && $CALLERBIN basecaller $CALLERPARAMS $MODELPARAMS inputs 2>> $ol 1> tmp.bam && samtools view -h tmp.bam|samtools fastq -n - | pigz 1> $oc 2>> $ol && rm -rf tmp.bam
    """
}

workflow BASECALL{ 
    take: collection

    main:

    SIGNALS = SAMPLES.collect{
        element -> return [element, BASECALLMAP[element].collect{ p -> return "${workflow.workDir}/../"+p }]
    }

    signals_ch = Channel.fromList(SIGNALS)
    
    dorado(signals_ch)

    emit:
    fastq = dorado.out.fastq
    logs = dorado.out.log
}