#!/usr/bin/env nextflow

//includes
//include {} from '../NextSnakes/lib/Collection.groovy'

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.preview.dsl=2

//Params from CL
REFERENCE = "${workflow.workDir}/../"+params.REFERENCE
REFDIR = "${workflow.workDir}/../"+params.REFDIR
BINS = "${workflow.workDir}/../"+params.BINS
UIDX = params.UIDX ?: null
INDEX = params.INDEX ?: null
INDEX2 = params.INDEX2 ?: null
PREFIX = params.PREFIX ?: null
THREADS = params.MAXTHREAD ?: null
PAIRED = params.PAIRED ?: null
RUNDEDUP = params.RUNDEDUP ?: null
STRANDED = params.STRANDED ?: null
IP = params.IP ?: null
CONDITION = params.CONDITION
COMBO = params.COMBO ?: ""
SCOMBO = params.SCOMBO ?: ""
ANNO = params.ANNOTATION ?: ""
SAMPLES = params.SAMPLES.split(',')
LONGSAMPLES = params.LONGSAMPLES.split(',')

//dummy
dummy = Channel.fromPath("${workflow.workDir}/../LOGS/NextSnakes.log")
