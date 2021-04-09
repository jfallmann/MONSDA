#!/usr/bin/env nextflow

//includes
//include {} from '../nextsnakes/lib/Collection.groovy'

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.preview.dsl=2

//Params from CL
REFERENCE = "${workflow.workDir}/../"+params.REFERENCE
REFDIR = "${workflow.workDir}/../"+params.REFDIR
UIDX = params.UIDX
INDEX = params.INDEX
INDEX2 = params.INDEX2
PREFIX = params.PREFIX
BINS = "${workflow.workDir}/../"+params.BINS
THREADS = params.MAXTHREAD
PAIRED = params.PAIRED
RUNDEDUP = params.RUNDEDUP
STRANDED = params.STRANDED
IP = params.IP
CONDITION = params.CONDITION
COMBO = params.COMBO
SCOMBO = params.SCOMBO

ANNO = params.ANNOTATION
SAMPLES = params.SAMPLES.split(',')
LONGSAMPLES = params.LONGSAMPLES.split(',')

//dummy
dummy = Channel.fromPath("${workflow.workDir}/../LOGS/RunNextflow.log")
