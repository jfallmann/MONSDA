#!/usr/bin/env nextflow

//includes
//include {} from '../nextsnakes/lib/Collection.groovy'

//Version Check
nextflowVersion = '>=20.01.0.5264'
nextflow.preview.dsl=2

//Params from CL
REFERENCE = params.REFERENCE
GENOME = params.GENOME
NAME = params.NAME
BINS = params.BINS
THREADS = params.MAXTHREAD
SOURCE = params.SOURCE
PAIRED = params.PAIRED
STRANDED = params.STRANDED
CONDITION = params.CONDITION
SAMPLES = params.SAMPLES.split(',')
LONGSAMPLES = params.LONGSAMPLES.split(',')

//dummy
dummy = null
