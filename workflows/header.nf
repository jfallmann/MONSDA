#!/usr/bin/env nextflow

//includes
//include {} from '../nextsnakes/lib/Collection.groovy'

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.preview.dsl=2

//Params from CL
REFERENCE = "${workflow.workDir}/../"+params.REFERENCE
REFDIR = "${workflow.workDir}/../"+params.REFDIR
GENOME = params.GENOME
INDEX = params.INDEX
PREFIX = params.PREFIX
BINS = "${workflow.workDir}/../"+params.BINS
THREADS = params.MAXTHREAD
PAIRED = params.PAIRED
DEDUP = params.DEDUP
STRANDED = params.STRANDED
IP = params.IP
CONDITION = params.CONDITION
SCONDITION = params.SCONDITION
SETS = params.SETS

ANNO = params.ANNOTATION
SAMPLES = params.SAMPLES.split(',')
LONGSAMPLES = params.LONGSAMPLES.split(',')

//dummy
dummy = Channel.fromPath("${workflow.workDir}/../LOGS/RunNextflow.log")
