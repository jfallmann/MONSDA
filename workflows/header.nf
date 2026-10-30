#!/usr/bin/env nextflow

// includes
// include {} from "../lib/Collection.groovy" 

// NOTE
// ALWAYS COMMENT LINES WITH '//', DO NOT USE MULTI LINE COMMENTS AS THE PARSER WILL NOT IGNORE MIDDLE LINES AND THIS WILL CAUSE CHAOS

//Version Check
nextflow.enable.dsl=2
//nextflowVersion = '>=20.01.0.5264'

//define unset Params
def get_always(parameter){
    return params.containsKey(parameter) ? params[parameter] : null
}

//Channel element order is arrival order, not read order, so R1/R2 of a collated pair must
//never be taken as [0]/[1] directly. Sorting by file name puts _R1 before _R2 deterministically.
def sort_reads(rds){
    def rdslist = (rds instanceof Collection) ? rds.toList() : [rds]
    return rdslist.toSorted{ a, b -> a.toString().split('/')[-1] <=> b.toString().split('/')[-1] }
}

def reads_for_mode(reads, mode){
    if (mode == 'paired' || !(PAIRED in ['paired', 'singlecell'])) return reads
    return reads.map{ item ->
        def files = (item instanceof Collection) ? item.toList() : [item]
        def read1 = files.find{ read ->
            read.toString().contains('_R1.') || read.toString().contains('_R1_')
        }
        if (read1 == null) throw new Exception("reads_for_mode: no R1 found in "+item)
        return read1
    }
}

def fastq_patterns(samples, mode){
    return samples.collect{ sample ->
        if (mode == 'paired') return "${workflow.workDir}/../FASTQ/"+sample+"_{R1,R2}.*fastq.gz"
        if (PAIRED in ['paired', 'singlecell']) return "${workflow.workDir}/../FASTQ/"+sample+"_R1.*fastq.gz"
        return "${workflow.workDir}/../FASTQ/"+sample+".*fastq.gz"
    }
}

//build_count_table.py zips the '-r' list positionally with '-c'/'-t'/'-b', so the replicate
//order must always follow the config derived REPS list and must never be taken from the
//order in which staged files arrive on a channel. This rewrites the '-r' list in place with
//the names of the staged files, keeping the configured order, and aborts on any mismatch.
def rep_basename(f){
    return f.toString().split('/')[-1]
}

def rep_key(f){
    return rep_basename(f).replaceFirst(/\.counts\.gz$/, '').replaceFirst(/_dedup$/, '')
}

def rep_args(repsargs, staged){
    def mt = (repsargs =~ /(^|\s)-r\s+(\S+)/)
    if (!mt.find()){
        throw new Exception("rep_args: no '-r' list found in REPS arguments '"+repsargs+"'")
    }
    def have = [:]
    def stagedlist = (staged instanceof Collection) ? staged.toList() : [staged]
    stagedlist.each{ f -> have[rep_key(f)] = rep_basename(f) }
    def ordered = mt.group(2).split(',').collect{ w ->
        def k = rep_key(w.trim())
        if (!have.containsKey(k)){
            throw new Exception("rep_args: no staged count file for replicate '"+w.trim()+"' (key '"+k+"'), staged keys: "+have.keySet())
        }
        return have[k]
    }
    if (ordered.size() != stagedlist.size()){
        throw new Exception("rep_args: "+stagedlist.size()+" staged count files but "+ordered.size()+" replicates in REPS arguments '"+repsargs+"'")
    }
    return repsargs.replaceFirst(/(^|\s)-r\s+\S+/, '$1-r '+ordered.join(','))
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
// Per-step pairedness/strandedness, default to the global values above and
// overridden by the workflow-specific MONSDA config key SEQUENCING, e.g. a
// quantseq MAPPING step can run 'single' while TRIMMING/DEDUP still see the
// global 'paired' mode needed for UMI extraction from read2.
QCPAIRED = get_always('QCPAIRED') ?: PAIRED
QCSTRANDED = params.containsKey('QCSTRANDED') ? get_always('QCSTRANDED') : STRANDED
TRIMMINGPAIRED = get_always('TRIMMINGPAIRED') ?: PAIRED
TRIMMINGSTRANDED = params.containsKey('TRIMMINGSTRANDED') ? get_always('TRIMMINGSTRANDED') : STRANDED
DEDUPPAIRED = get_always('DEDUPPAIRED') ?: PAIRED
DEDUPSTRANDED = params.containsKey('DEDUPSTRANDED') ? get_always('DEDUPSTRANDED') : STRANDED
MAPPINGPAIRED = get_always('MAPPINGPAIRED') ?: PAIRED
MAPPINGSTRANDED = get_always('MAPPINGSTRANDED') ?: STRANDED
PEAKPAIRED = get_always('PEAKPAIRED') ?: MAPPINGPAIRED
PEAKSTRANDED = get_always('PEAKSTRANDED') ?: MAPPINGSTRANDED
TRACKSPAIRED = get_always('TRACKSPAIRED') ?: MAPPINGPAIRED
TRACKSSTRANDED = get_always('TRACKSSTRANDED') ?: MAPPINGSTRANDED
COUNTINGPAIRED = get_always('COUNTINGPAIRED') ?: MAPPINGPAIRED
COUNTINGSTRANDED = get_always('COUNTINGSTRANDED') ?: MAPPINGSTRANDED
CIRCSPAIRED = get_always('CIRCSPAIRED') ?: MAPPINGPAIRED
CIRCSSTRANDED = get_always('CIRCSSTRANDED') ?: MAPPINGSTRANDED
FUSIONSPAIRED = get_always('FUSIONSPAIRED') ?: MAPPINGPAIRED
FUSIONSSTRANDED = get_always('FUSIONSSTRANDED') ?: MAPPINGSTRANDED
DEPAIRED = get_always('DEPAIRED') ?: MAPPINGPAIRED
DESTRANDED = get_always('DESTRANDED') ?: MAPPINGSTRANDED
DEUPAIRED = get_always('DEUPAIRED') ?: MAPPINGPAIRED
DEUSTRANDED = get_always('DEUSTRANDED') ?: MAPPINGSTRANDED
DASPAIRED = get_always('DASPAIRED') ?: MAPPINGPAIRED
DASSTRANDED = get_always('DASSTRANDED') ?: MAPPINGSTRANDED
// DTU's default salmon path quantifies directly from TRIMMED_FASTQ, independent
// of the MAPPING step, so it falls back to the global layout, not MAPPING's
DTUPAIRED = get_always('DTUPAIRED') ?: PAIRED
DTUSTRANDED = get_always('DTUSTRANDED') ?: STRANDED
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
RSAMPLES = fastq_patterns(SAMPLES, PAIRED == 'singlecell' ? 'paired' : PAIRED)
samples_ch = Channel.fromPath(RSAMPLES)