#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def get_always(parameter){
    return params.containsKey(parameter) ? params[parameter] : null
}
params.gREFERENCE = "${workflow.workDir}/../"+get_always('REFERENCE')
params.gREFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
params.gBINS = get_always('BINS')
params.gTHREADS = get_always('MAXTHREAD')
params.gPAIRED = get_always('PAIRED') ?: null
params.gRUNDEDUP = get_always('RUNDEDUP') ?: null
params.gPREDEDUP = get_always('PREDEDUP') ?: null
params.gSTRANDED = get_always('STRANDED') ?: null
params.gIP = get_always('IP') ?: null
params.gCONDITION = get_always('CONDITION') ?: null
params.gCOMBO = get_always('COMBO') ?: ''
params.gSCOMBO = get_always('SCOMBO') ?: ''
params.gSAMPLES = get_always('SAMPLES').split(',') ?: null
params.gLONGSAMPLES = get_always('LONGSAMPLES').split(',') ?: null
params.gSHORTSAMPLES = get_always('SHORTSAMPLES').split(',') ?: null
params.gSETS = get_always('SETS') ?: null
def dummy(){
    def dummy = Channel.fromPath("${workflow.workDir}/../LOGS/MONSDA.log")
    return dummy
}
def samples_ch(){
    def samples_ch = Channel.fromPath(params.gRSAMPLES)
    return samples_ch
}
params.gQCENV = get_always('PREQCENV')
params.gQCBIN = get_always('PREQCBIN')
params.gQCPARAMS = get_always('fastqc_params_QC') ?: ''
params.gTRIMENV = get_always('TRIMMINGENV')
params.gTRIMBIN = get_always('TRIMMINGBIN')
params.gTRIMPARAMS = get_always('cutadapt_params_TRIM') ?: ''
params.gDEDUPENV = get_always('DEDUPENV')
params.gDEDUPBIN = get_always('DEDUPBIN')
params.gDEDUPPARAMS = get_always('picard_params_DEDUP') ?: ''
params.gJAVAPARAMS = get_always('picard_params_JAVA') ?: ''
params.gMAPENV = get_always('MAPPINGENV')
params.gMAPBIN = get_always('MAPPINGBIN')
params.gMAPIDX = get_always('MAPPINGIDX')
params.gMAPUIDX = { def MAPUIDX = get_always('MAPPINGUIDX'); MAPUIDX = MAPUIDX.replace('.idx',''); return MAPUIDX }()
params.gMAPUIDXNAME = get_always('MAPPINGUIDXNAME')
params.gMAPREF = get_always('MAPPINGREF')
params.gMAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
params.gMAPANNO = get_always('MAPPINGANNO')
params.gMAPPREFIX = get_always('MAPPINGPREFIX')
params.gIDXPARAMS = get_always('star_params_INDEX') ?: ''
params.gMAPPARAMS = get_always('star_params_MAP') ?: ''
params.gMQCENV = get_always('PREQCENV')
params.gMQCBIN = get_always('PREQCBIN')
params.gMQCPARAMS = get_always('fastqc_params_MULTI') ?: ''
params.gRSAMPLES = {
if (params.gPAIRED == 'paired' || params.gPAIRED == 'singlecell'){
    return params.gSAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+"_{R2,R1}.*fastq.gz"
    }
}else{
    return params.gSAMPLES.collect{
        element -> return "${workflow.workDir}/../FASTQ/"+element+".*fastq.gz"
    }
}
}()

process qc_raw{
    conda "<REPO>/envs/${params.gQCENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gQCENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} ${params.gQCPARAMS} --noextract -f fastq $read
    """
}
workflow QC_RAW{
    take:
    collection
    main:
    qc_raw(samples_ch())
    emit:
    qc = qc_raw.out.fastqc_results
}
process qc_trimmed{
    conda "<REPO>/envs/${params.gQCENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gQCENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} ${params.gQCPARAMS} --noextract -f fastq $read
    """
}
workflow QC_TRIMMING{
    take: collection
    main:
    qc_trimmed(collection)
    emit:
    qc = qc_trimmed.out.fastqc_results
}
process qc_mapped{
    conda "<REPO>/envs/${params.gQCENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gQCENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path map
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} ${params.gQCPARAMS} -f bam $map
    """
}
workflow QC_MAPPING{
    take: collection
    main:
    qc_mapped(collection)
    emit:
    qc = qc_mapped.out.fastqc_results
}
process qc_dedup{
    conda "<REPO>/envs/${params.gQCENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gQCENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} ${params.gQCPARAMS} --noextract -f fastq $read
    """
}
workflow QC_DEDUP{
    take: collection
    main:
    qc_dedup(collection)
    emit:
    qc = qc_dedup.out.fastqc_results
}


process trim{
    conda "<REPO>/envs/${params.gTRIMENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gTRIMENV}"+"-VERSION"
    cpus 4
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_trimmed.fastq.gz") > 0)                "TRIMMED_FASTQ/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed|_dedup/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}"
        else if (filename.indexOf(".log") >0)              "LOGS/${params.gCOMBO}/${params.gCONDITION}/TRIMMING/cutadapt/${file(filename).getSimpleName()}.log"
        else null
    }
    input:
    path reads
    output:
    path "*_trimmed.fastq.gz", emit: trim
    path "*trimming_report.txt", emit: rep
    script:
    if (params.gPAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        o = file(r1).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"
        p = file(r2).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"
        r = file(r1).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimming_report.txt"
        """
        ${params.gTRIMBIN} --cores ${task.cpus} ${params.gTRIMPARAMS} -o $o -p $p $r1 $r2 &> $r
        """
    }
    else{
        o = file(reads).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"
        r = file(reads).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimming_report.txt"
        """
        ${params.gTRIMBIN} --cores ${task.cpus} ${params.gTRIMPARAMS} -o $o $reads &> $r
        """
    }
}
workflow TRIMMING{
    take: 
    collection    
    main:
    if ( params.gPREDEDUP == 'enabled' ){  
        trim(collection)
    }else {        
        if (params.gPAIRED == 'paired'){
            trim(samples_ch().collate(2))
        } else{
            trim(samples_ch().collate(1))
        }
    }
    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}


process dedup_bam{
    conda "<REPO>/envs/${params.gDEDUPENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gDEDUPENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))              "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_dedup.bam.bai") > 0)  "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("dedup.log") > 0)       "LOGS/${params.gCOMBO}/${params.gCONDITION}/DEDUP/picard/${file(filename).getName()}"
        else if (filename.indexOf("metrix.txt") > 0)      "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
        else null
    }
    publishDir "${workflow.workDir}/../LOGS/${params.gCOMBO}/${params.gCONDITION}/DEDUP/picard" , mode: 'copy', pattern: "*_dedup_metrix.txt"
    input:
    path todedup
    path bami
    output:
    path "*_dedup.bam", emit: bam
    path "*_dedup.bam.bai", emit: bai
    path "*_dedup.log", emit: logs
    path "*_dedup_metrix.txt", emit: metrics
    script:
    bams = todedup[0]
    bais = todedup[1]
    outf = bams.getSimpleName()+"_dedup.bam"
    outl = bams.getSimpleName()+"_dedup.log"
    outm = bams.getSimpleName()+"_dedup_metrix.txt"
    """
    mkdir -p TMP && ${params.gDEDUPBIN} ${params.gJAVAPARAMS} MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR TMP --INPUT $bams --OUTPUT $outf --METRICS_FILE $outm ${params.gDEDUPPARAMS} &> $outl && samtools index $outf &>> $outl
    """
}
workflow DEDUPBAM{
    take:
    map
    mapi
    mapu
    mapui
    main:   
    dedup_bam(map.concat(mapu), mapi.concat(mapui))
    emit:
    dedup = dedup_bam.out.bam
    dedupbai = dedup_bam.out.bai
    deduplog = dedup_bam.out.logs
}


process collect_tomap{
    input:
    path check
    output:
    path "collect.txt", emit: done
    script:
    """
    echo "$check Collection successful!" > collect.txt
    """
}
process star_idx{
    conda "<REPO>/envs/${params.gMAPENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gMAPENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    label 'big_mem'
    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)             "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/star/index.log"
        else if (filename.indexOf(".idx") > 0)           "${params.gMAPIDX}"
        else                                             "${params.gMAPUIDX}"
    }
    input:
    path genome
    path anno
    output:
    path "${params.gMAPUIDXNAME}", emit: idx
    path "*.out", emit: idxlog
    path "*.idx", emit: tmpidx
    script:
    gen =  genome.getName()
    an  = anno.getName()
    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p ${params.gMAPUIDXNAME} && rm -rf STARTMP && ${params.gMAPBIN} ${params.gIDXPARAMS} --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir ${params.gMAPUIDXNAME} --genomeFastaFiles tmp.fa && touch ${params.gMAPUIDXNAME} && ln -s ${params.gMAPUIDXNAME} star.idx && rm -f tmp.fa tmp_anno && ln -fs ${params.gMAPUIDXNAME}/* . && rm -rf STARTMP
    """
}
process star_mapping{
    conda "<REPO>/envs/${params.gMAPENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gMAPENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    label 'big_mem'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped") > 0)       "UNMAPPED/${params.gCOMBO}/${params.gCONDITION}/"+"${file(filename).getName()}"
        else if (filename.indexOf(".log") >0)        "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/star/"+"${file(filename).getName()}"
        else if (filename.indexOf(".out") >0)        "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/star/"+"${file(filename).getName()}"
        else if (filename.indexOf(".tab") >0)        "MAPPED/${params.gCOMBO}/${params.gCONDITION}/"+"${filename}"
        else                                         "MAPPED/${params.gCOMBO}/${params.gCONDITION}/"+"${filename}"
    }
    input:
    path reads
    output:
    path "*_mapped.sam.gz", emit: maps
    path "*.out", emit: outs
    path "*.log", emit: logs
    path "*Log*.out", emit: xtralogs
    path "*.tab", emit: sjtab
    path "*_unmapped.fastq.gz", includeInputs:false, emit: unmapped
    script:
    idx = reads[0]
    idxdir = idx.toRealPath()
    if (params.gPAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        a = "Trimming_report.txt"
        fn = file(r1).getSimpleName().replaceAll(/_R1(_dedup)?_trimmed$/,"")
        of = fn+'.Aligned.out.sam'
        gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
        """
        ${params.gMAPBIN} ${params.gMAPPARAMS} --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2 --outFileNamePrefix ${fn}. --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && touch ${fn}.Unmapped.out.mate1 ${fn}.Unmapped.out.mate2 && cat ${fn}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R1_unmapped.fastq.gz && cat ${fn}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R2_unmapped.fastq.gz && mv *Log.out ${fn}_mapping.log
        """
    }
    else{
        if (params.gPAIRED != 'singlecell'){
            read = reads[1]
            fn = file(reads[1]).getSimpleName().replaceAll(/(_dedup)?_trimmed$/,"")+"."
            of = fn+'Aligned.out.sam'
            gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
            """
            ${params.gMAPBIN} ${params.gMAPPARAMS} --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $read --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed -r 's/\\.Unmapped.out.mate1.gz/_unmapped.fastq.gz/')"; done && mv *Log.out ${fn}_mapping.log
            """
        }
        else{
            if (params.gSTRANDED == 'fr'){
                stranded = '--soloStrand Forward'
            }else if (params.gSTRANDED == 'rf'){
                stranded = '--soloStrand Reverse'
            }else{
                stranded = '--soloStrand Unstranded'
            }
            r1 = reads[1]
            fn = file(r1).getSimpleName().replaceAll(/_R1(_dedup)?_trimmed$/,"")
            r2 = "${workflow.workDir}/../FASTQ/${params.gCONDITION}/"+file(reads[2]).getSimpleName().replaceAll(/\QR2_trimmed\E/,"R2.fastq.gz")
            if (params.gMAPPARAMS.contains('--soloBarcodeMate 1')){
                t = r2
                r2 = r1
                r1 = t
            }
            of = fn+'.Aligned.sortedByCoord.out.bam'
            gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")
            """
            ${params.gMAPBIN} ${params.gMAPPARAMS} $stranded --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2  --outFileNamePrefix ${fn}. --outReadsUnmapped Fastx && samtools view -h ${of} | gzip > $gf && rm -f $of && touch ${fn}.Unmapped.out.mate1 ${fn}.Unmapped.out.mate2 && cat ${fn}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R1_unmapped.fastq.gz && at ${fn}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R2_unmapped.fastq.gz && mv *Log.out ${fn}_mapping.log
            """
        }
    }
}
workflow MAPPING{
    take: collection
    main:
    checkidx = file(params.gMAPUIDX)
    if (checkidx.exists()){
        idxfile = Channel.fromPath(params.gMAPUIDX)
        star_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(params.gMAPREF)
        annofile = Channel.fromPath(params.gMAPANNO)
        star_idx(genomefile, annofile)
        star_mapping(star_idx.out.idx.combine(collection))
    }
    emit:
    mapped  = star_mapping.out.maps
    logs = star_mapping.out.logs
}


process sortsam{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"samtools"+"-VERSION"
    cpus params.gTHREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sam.gz") > 0)     "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.sam.gz"
        else null
    }
    input:
    path map
    output:
    path "*_mapped_sorted.sam.gz", includeInputs:false, emit: sam
    script:
    fn = file(map[0]).getSimpleName()
    sorted = fn+'_sorted.sam.gz'
    def sortmem = Math.ceil(task.memory.giga as double) as int 
    """
    set +o pipefail; mkdir -p TMP; samtools view -H $map|grep -P '^@HD' |pigz -p ${task.cpus} -f > tmphead; samtools view -H $map|grep -P '^@SQ'|sort -t\$'\\t' -k1,1 -k2,2V |pigz -p ${task.cpus} -f >> tmphead ; samtools view -H $map|grep -P '^@RG'|pigz -p ${task.cpus} -f >> tmphead ; samtools view -H $map|grep -P '^@PG'|pigz -p ${task.cpus} -f >> tmphead ; export LC_ALL=C;zcat $map | grep -v \"^@\"|sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k3,3V -k4,4n - |pigz -p ${task.cpus} -f > tmpfile; cat tmphead tmpfile > $sorted && rm -f tmphead tmpfile
    """
}
process sam2bam{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"samtools"+"-VERSION"
    cpus params.gTHREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/samtools/${file(filename).getSimpleName()}.log"
        else null
    }
    input:
    path sam
    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: logs
    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+".bam"
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads ${task.cpus} && samtools index $bam 2> sam2bam.log
    """
}
process uniqsam{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"samtools"+"-VERSION"
    cpus params.gTHREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("unique.sam.gz") > 0)   "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.sam.gz"
        else if (filename.indexOf(".log") > 0)       "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/samtools/${file(filename).getSimpleName()}.log"
        else null
    }
    input:
    path sam
    output:
    path "*_unique.sam.gz", includeInputs:false, emit: sam
    path "*.log", emit: logs
    script:
    fn = file(sam[0]).getSimpleName()
    uniq = fn+'_unique.sam.gz'
    if (params.gCOMBO.indexOf('-bwa') > 0){
        bwa = 'bwa'
    } else{
        bwa = ''
    }
    """
    ${params.gBINS}/Shells/UniqueSam_woPicard.sh $sam $uniq ${task.cpus} $bwa 2> uniq.log
    """
}
process sam2bamuniq{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"samtools"+"-VERSION"
    cpus params.gTHREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${params.gCOMBO}/${params.gCONDITION}/MAPPING/samtools/${file(filename).getSimpleName()}.log"
        else null
    }
    input:
    path sam
    output:
    path "*.bam", emit: bam
    path "*.bai", emit: bai
    path "*.log", emit: logs
    script:
    fn = file(sam[0]).getSimpleName()
    bam = fn+'.bam'
    """
    zcat $sam | samtools view -bS - | samtools sort -T $fn -o $bam --threads ${task.cpus} && samtools index $bam 2> uniquebam.log
    """
}
workflow POSTMAPPING{
    take: collection
    main:
    sortsam(collection)
    sam2bam(sortsam.out.sam)
    uniqsam(sortsam.out.sam)
    sam2bamuniq(uniqsam.out.sam)
    emit:
    postmap  = sam2bam.out.bam
    postbai  = sam2bam.out.bai
    postmapuni = sam2bamuniq.out.bam
    postunibai = sam2bamuniq.out.bai
}


process mqc{
    conda "<REPO>/envs/${params.gMQCENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gMQCENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getSimpleName()}.html"
        else "QC/Multi/${params.gCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
    }
    input:
    path others, stageAs: 'mqc_input??/*'
    output:
    path "*.zip", emit: mqc
    path "*.html", emit: html
    script:
    """
    touch $others
    OUT=\${PWD}
    LOGDIR="${workflow.workDir}/../LOGS"
    SCAN="\${LOGDIR}/${params.gCOMBO}/${params.gCONDITION}"
    COLLECT="\${SCAN}/MULTIQC/collect"
    VERSIONS="\${LOGDIR}/versions.txt"
    EXTRA=""
    BASE_QC_DIR="${workflow.workDir}/../QC"
    COMBO_VAL="${params.gCOMBO}"
    CONDITION_VAL="${params.gCONDITION}"
    mkdir -p "\$COLLECT"
    for i in $others; do
        cp -f "\$i" "\$COLLECT"/
    done
    # If the corresponding fastqc combo exists, include its output in the MultiQC report.
    FQ_COMBO="\${COMBO_VAL/rustqc/fastqc}"
    FQ_DIR="\${BASE_QC_DIR}/\${FQ_COMBO}/\${CONDITION_VAL}"
    if [[ "\$FQ_COMBO" != "\$COMBO_VAL" && -d "\$FQ_DIR" ]]; then
        EXTRA="\$FQ_DIR"
    fi
    MODS=""
    if [[ -f "\$VERSIONS" ]]; then
        MODS=\$(grep -v '^#' "\$VERSIONS" | cut -f3 | tr ',' '\\n' | grep -vx '-' | sort -u | sed 's/^/-m /' | tr '\\n' ' ')
        cp -f "\$VERSIONS" "\$OUT"/
    fi
    export LC_ALL=C.UTF-8
    multiqc -f ${params.gMQCPARAMS} \$MODS -k json -z -s -o "\$OUT" "\$SCAN" \$EXTRA
    """
}
workflow MULTIQC{
    take:
    otherqcs
    main:
    mqc(otherqcs.collect())
    emit:
    mqcres = mqc.out.mqc
}




workflow {
    QC_RAW(dummy())
    TRIMMING(dummy())
    QC_TRIMMING(TRIMMING.out.trimmed)
    MAPPING(TRIMMING.out.trimmed)
    POSTMAPPING(MAPPING.out.mapped)
    DEDUPBAM(POSTMAPPING.out.postmap, POSTMAPPING.out.postbai, POSTMAPPING.out.postmapuni, POSTMAPPING.out.postunibai)
    QC_MAPPING(POSTMAPPING.out.postmap.concat(POSTMAPPING.out.postmapuni.concat(DEDUPBAM.out.dedup)))
    MULTIQC(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc.concat(QC_MAPPING.out.qc.concat(MAPPING.out.logs))).collect())
}

