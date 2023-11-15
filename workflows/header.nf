#!/usr/bin/env nextflow

// includes
// include {} from "../lib/Collection.groovy" 

// NOTE
// ALWAYS COMMENT LINES WITH '//', DO NOT USE MULTI LINE COMMENTS AS THE PARSER WILL NOT IGNORE MIDDLE LINES AND THIS WILL CAUSE CHAOS

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.enable.dsl=2

//define unset Params
def get_always(parameter){
    if (!params.containsKey(parameter)){
        params.put(parameter, null)
    }
    return params[parameter]
}

//Params from CL
REFERENCE = "${workflow.workDir}/../"+get_always('REFERENCE')
REFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
BINS = get_always('BINS')
THREADS = get_always('MAXTHREAD')
PAIRED = get_always('PAIRED') ?: null
RUNDEDUP = get_always('RUNDEDUP') ?: null
PREDEDUP = get_always('PREDEDUP') ?: null
STRANDED = get_always('STRANDED') ?: null
IP = get_always('IP') ?: null
CONDITION = get_always('CONDITION') ?: null
COMBO = get_always('COMBO') ?: ''
SCOMBO = get_always('SCOMBO') ?: ''
SAMPLES = get_always('SAMPLES').split(',') ?: null
LONGSAMPLES = get_always('LONGSAMPLES').split(',') ?: null
SHORTSAMPLES = get_always('SHORTSAMPLES').split(',') ?: null
SETS = get_always('SETS') ?: null
//dummy
dummy = Channel.fromPath("${workflow.workDir}/../LOGS/MONSDA.log")

//SAMPLE CHANNELS
if (PAIRED == 'paired' || PAIRED == 'singlecell'){
    RSAMPLES = SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
    }
}else{
    RSAMPLES=SAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
    }
}

samples_ch = Channel.fromPath(RSAMPLES)