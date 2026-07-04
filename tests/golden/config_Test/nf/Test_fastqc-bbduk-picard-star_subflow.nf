#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def get_always(parameter){
    if (!params.containsKey(parameter)){
        params.put(parameter, null)
    }
    return params[parameter]
}
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
dummy = Channel.fromPath("${workflow.workDir}/../LOGS/MONSDA.log")
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

QCENV=get_always('PREQCENV')
QCBIN=get_always('PREQCBIN')
QCPARAMS = get_always('fastqc_params_QC') ?: ''
process qc_raw{
    conda "<REPO>/envs/$QCENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} $QCPARAMS --noextract -f fastq $read
    """
}
workflow QC_RAW{
    take:
    collection
    main:
    qc_raw(samples_ch)
    emit:
    qc = qc_raw.out.fastqc_results
}
process qc_trimmed{
    conda "<REPO>/envs/$QCENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} $QCPARAMS --noextract -f fastq $read
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
    conda "<REPO>/envs/$QCENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path map
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} $QCPARAMS -f bam $map
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
    conda "<REPO>/envs/$QCENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else null
    }
    input:
    path read
    output:
    path "*.{zip,html}", emit: fastqc_results
    script:
    """
    fastqc --quiet -t ${task.cpus} $QCPARAMS --noextract -f fastq $read
    """
}
workflow QC_DEDUP{
    take: collection
    main:
    qc_dedup(collection)
    emit:
    qc = qc_dedup.out.fastqc_results
}


TRIMENV=get_always('TRIMMINGENV')
TRIMBIN=get_always('TRIMMINGBIN')
TRIMPARAMS = get_always('bbduk_params_TRIM') ?: ''
process trim{
    conda "<REPO>/envs/$TRIMENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$TRIMENV"
    cpus 4
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_trimmed.fastq.gz") > 0)     "TRIMMED_FASTQ/${COMBO}/${CONDITION}/${file(filename).getSimpleName().replaceAll(/_val_\d{1}|_trimmed|_dedup/,"")}_trimmed.fastq.gz"
        else if (filename.indexOf("report.txt") >0)        "TRIMMED_FASTQ/${COMBO}/${CONDITION}/${file(filename).getSimpleName().replaceAll(/.fastq.gz/,"")}_trimming_report.txt"
        else if (filename.indexOf(".log") >0)              "LOGS/${COMBO}/${CONDITION}/TRIMMING/${file(filename).getSimpleName()}.log"
        else null
    }
    input:
    path reads
    output:
    path "*_trimmed.fastq.gz", emit: trim
    path "*trimming_report.txt", emit: rep
    script:
    if (PAIRED == 'paired'){
        r1 = reads[0]
        r2 = reads[1]
        o = file(r1).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"
        p = file(r2).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"
        r = file(r1).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimming_report.txt"
        """
        $TRIMBIN $TRIMPARAMS t=${task.cpus} in1=$r1 in2=$r2 out1=$o out2=$p &> $r 
        """
    }
    else{
        o = file(reads).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimmed.fastq.gz"        
        r = file(reads).getSimpleName().replaceAll(/_dedup/,"").replaceAll(/.fastq.gz/,"")+"_trimming_report.txt"
        """
        $TRIMBIN $TRIMPARAMS t=${task.cpus} in=$r1 out=$o &> $r 
        """
    }
}
workflow TRIMMING{
    take: 
    collection    
    main:
    if ( PREDEDUP == 'enabled' ){  
        trim(collection)
    }else {        
        if (PAIRED == 'paired'){
            trim(samples_ch.collate(2))
        } else{
            trim(samples_ch.collate(1))
        }
    }
    emit:
    trimmed = trim.out.trim
    report  = trim.out.rep
}


DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')
DEDUPPARAMS = get_always('picard_params_DEDUP') ?: ''
JAVAPARAMS = get_always('picard_params_JAVA') ?: ''
process dedup_bam{
    conda "<REPO>/envs/$DEDUPENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$DEDUPENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))              "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_dedup.bam.bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("dedup.log") > 0)       "LOGS/${COMBO}/${CONDITION}/DEDUP/${file(filename).getName()}"
        else if (filename.indexOf("metrix.txt") > 0)      "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else null
    }
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
    mkdir -p TMP && $DEDUPBIN $JAVAPARAMS MarkDuplicates --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER coordinate --TMP_DIR TMP --INPUT $bams --OUTPUT $outf --METRICS_FILE $outm $DEDUPPARAMS &> $outl && samtools index $outf &>> $outl
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


MAPENV = get_always('MAPPINGENV')
MAPBIN = get_always('MAPPINGBIN')
MAPIDX = get_always('MAPPINGIDX')
MAPUIDX = get_always('MAPPINGUIDX')
MAPUIDXNAME = get_always('MAPPINGUIDXNAME')
MAPREF = get_always('MAPPINGREF')
MAPREFDIR = "${workflow.workDir}/../"+get_always('MAPPINGREFDIR')
MAPANNO = get_always('MAPPINGANNO')
MAPPREFIX = get_always('MAPPINGPREFIX')
MAPUIDX = MAPUIDX.replace('.idx','')
IDXPARAMS = get_always('star_params_INDEX') ?: ''
MAPPARAMS = get_always('star_params_MAP') ?: ''
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
    conda "<REPO>/envs/$MAPENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow', overwrite: true,
    saveAs: {filename ->
        if (filename.indexOf("Log.out") > 0)             "LOGS/${COMBO}/${CONDITION}/star_index.log"
        else if (filename.indexOf(".idx") > 0)           "$MAPIDX"
        else                                             "$MAPUIDX"
    }
    input:
    path genome
    path anno
    output:
    path "$MAPUIDXNAME", emit: idx
    path "*.out", emit: idxlog
    path "*.idx", emit: tmpidx
    script:
    gen =  genome.getName()
    an  = anno.getName()
    """
    zcat $gen > tmp.fa && zcat $an > tmp_anno && mkdir -p $MAPUIDXNAME && $MAPBIN $IDXPARAMS --runThreadN ${task.cpus} --runMode genomeGenerate --outTmpDir STARTMP --genomeDir $MAPUIDXNAME --genomeFastaFiles tmp.fa && touch $MAPUIDXNAME && ln -s $MAPUIDXNAME star.idx && rm -f tmp.fa tmp_anno && ln -fs $MAPUIDXNAME/* . && rm -rf STARTMP
    """
}
process star_mapping{
    conda "<REPO>/envs/$MAPENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$MAPENV"
    cpus THREADS
	cache 'lenient'
    label 'big_mem'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_unmapped") > 0)       "UNMAPPED/${COMBO}/${CONDITION}/"+"${file(filename).getName()}"
        else if (filename.indexOf(".log") >0)        "LOGS/${COMBO}/${CONDITION}/MAPPING/star_"+"${file(filename).getName()}"
        else if (filename.indexOf(".out") >0)        "LOGS/${COMBO}/${CONDITION}/MAPPING/star_"+"${file(filename).getName()}"
        else if (filename.indexOf(".tab") >0)        "MAPPED/${COMBO}/${CONDITION}/"+"${filename}"
        else                                         "MAPPED/${COMBO}/${CONDITION}/"+"${filename}"
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
    if (PAIRED == 'paired'){
        r1 = reads[1]
        r2 = reads[2]
        a = "Trimming_report.txt"
        fn = file(r1).getSimpleName().replaceAll(/_R1(_dedup)?_trimmed$/,"")
        of = fn+'.Aligned.out.sam'
        gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
        """
        $MAPBIN $MAPPARAMS --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2 --outFileNamePrefix ${fn}. --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && touch ${fn}.Unmapped.out.mate1 ${fn}.Unmapped.out.mate2 && cat ${fn}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R1_unmapped.fastq.gz && cat ${fn}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R2_unmapped.fastq.gz && mv *Log.out ${fn}_mapping.log
        """
    }
    else{
        if (PAIRED != 'singlecell'){
            read = reads[1]
            fn = file(reads[1]).getSimpleName().replaceAll(/(_dedup)?_trimmed$/,"")+"."
            of = fn+'Aligned.out.sam'
            gf = of.replaceAll(/\Q.Aligned.out.sam\E/,"_mapped.sam.gz")
            """
            $MAPBIN $MAPPARAMS --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $read --outFileNamePrefix $fn --outReadsUnmapped Fastx && gzip -c $of > $gf && rm -f $of && gzip *Unmapped.out* && for f in *mate*.gz; do mv "\$f" "\$(echo "\$f" | sed -r 's/\\.Unmapped.out.mate1.gz/_unmapped.fastq.gz/')"; done && mv *Log.out ${fn}_mapping.log
            """
        }
        else{
            if (STRANDED == 'fr'){
                stranded = '--soloStrand Forward'
            }else if (STRANDED == 'rf'){
                stranded = '--soloStrand Reverse'
            }else{
                stranded = '--soloStrand Unstranded'
            }
            r1 = reads[1]
            fn = file(r1).getSimpleName().replaceAll(/_R1(_dedup)?_trimmed$/,"")
            r2 = "${workflow.workDir}/../FASTQ/${CONDITION}/"+file(reads[2]).getSimpleName().replaceAll(/\QR2_trimmed\E/,"R2.fastq.gz")
            if (MAPPARAMS.contains('--soloBarcodeMate 1')){
                t = r2
                r2 = r1
                r1 = t
            }
            of = fn+'.Aligned.sortedByCoord.out.bam'
            gf = of.replaceAll(/\Q.Aligned.sortedByCoord.out.bam\E/,"_mapped.sam.gz")
            """
            $MAPBIN $MAPPARAMS $stranded --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --runThreadN ${task.cpus} --genomeDir $idxdir --readFilesCommand zcat --readFilesIn $r1 $r2  --outFileNamePrefix ${fn}. --outReadsUnmapped Fastx && samtools view -h ${of} | gzip > $gf && rm -f $of && touch ${fn}.Unmapped.out.mate1 ${fn}.Unmapped.out.mate2 && cat ${fn}.Unmapped.out.mate1 | paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R1_unmapped.fastq.gz && at ${fn}.Unmapped.out.mate2| paste - - - - |tr \"\\t\" \"\\n\"| gzip > ${fn}_R2_unmapped.fastq.gz && mv *Log.out ${fn}_mapping.log
            """
        }
    }
}
workflow MAPPING{
    take: collection
    main:
    checkidx = file(MAPUIDX)
    if (checkidx.exists()){
        idxfile = Channel.fromPath(MAPUIDX)
        star_mapping(idxfile.combine(collection))
    }
    else{
        genomefile = Channel.fromPath(MAPREF)
        annofile = Channel.fromPath(MAPANNO)
        star_idx(genomefile, annofile)
        star_mapping(star_idx.out.idx.combine(collection))
    }
    emit:
    mapped  = star_mapping.out.maps
    logs = star_mapping.out.logs
}


process sortsam{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://docker.io/jfallmann/monsda:"+"samtools"
    cpus THREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sam.gz") > 0)     "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.sam.gz"
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
    container "oras://docker.io/jfallmann/monsda:"+"samtools"
    cpus THREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
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
    container "oras://docker.io/jfallmann/monsda:"+"samtools"
    cpus THREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("unique.sam.gz") > 0)   "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.sam.gz"
        else if (filename.indexOf(".log") > 0)       "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
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
    if (COMBO.indexOf('-bwa') > 0){
        bwa = 'bwa'
    } else{
        bwa = ''
    }
    """
    $BINS/Shells/UniqueSam_woPicard.sh $sam $uniq ${task.cpus} $bwa 2> uniq.log
    """
}
process sam2bamuniq{
    conda "<REPO>/envs/samtools.yaml"
    container "oras://docker.io/jfallmann/monsda:"+"samtools"
    cpus THREADS
    memory { 16.GB * (1 << ((task.attempt ?: 1) - 1)) }
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith(".bam"))       "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam"
        else if (filename.indexOf(".bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.bam.bai"
        else if (filename.indexOf(".log") > 0)  "LOGS/${COMBO}/${CONDITION}/MAPPING/${file(filename).getSimpleName()}.log"
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


QCENV=get_always('PREQCENV')
QCBIN=get_always('PREQCBIN')
QCPARAMS = get_always('fastqc_params_MULTI') ?: ''
process mqc{
    conda "<REPO>/envs/$QCENV"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"$QCENV"
    cpus THREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("zip") > 0)          "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.zip"
        else if (filename.indexOf("html") > 0)    "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getSimpleName()}.html"
        else "QC/Multi/${COMBO}/${CONDITION}/${file(filename).getName()}"
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
    LIST=multiqc_inputs.txt
    TMP_LIST=multiqc_inputs_unique.txt
    BASE_QC_DIR="${workflow.workDir}/../QC"
    COMBO_VAL="${COMBO}"
    CONDITION_VAL="${CONDITION}"
    for i in $others; do
        dirname "\$i" >> "\$LIST"
    done
    # If this is a rustqc combo and the corresponding fastqc combo exists,
    # include the fastqc output directory in the same MultiQC report.
    if [[ "\$COMBO_VAL" == *"rustqc"* ]]; then
        FQ_COMBO="\${COMBO_VAL/rustqc/fastqc}"
        FQ_DIR="\${BASE_QC_DIR}/\${FQ_COMBO}/\${CONDITION_VAL}"
        if [[ -d "\$FQ_DIR" ]]; then
            echo "\$FQ_DIR" >> "\$LIST"
        fi
    fi
    sort -u "\$LIST" > "\$TMP_LIST"
    export LC_ALL=en_US.utf8
    export LC_ALL=C.UTF-8
    multiqc -f --exclude picard --exclude gatk -k json -z -s -o "\$OUT" -l "\$TMP_LIST"
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
    QC_RAW(dummy)
    TRIMMING(dummy)
    QC_TRIMMING(TRIMMING.out.trimmed)
    MAPPING(TRIMMING.out.trimmed)
    POSTMAPPING(MAPPING.out.mapped)
    DEDUPBAM(POSTMAPPING.out.postmap, POSTMAPPING.out.postbai, POSTMAPPING.out.postmapuni, POSTMAPPING.out.postunibai)
    QC_MAPPING(POSTMAPPING.out.postmap.concat(POSTMAPPING.out.postmapuni.concat(DEDUPBAM.out.dedup)))
    MULTIQC(QC_RAW.out.qc.concat(QC_TRIMMING.out.qc.concat(QC_MAPPING.out.qc.concat(MAPPING.out.logs))).collect())
}

