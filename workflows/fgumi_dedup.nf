DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

DEDUPPARAMS = get_always('fgumi_params_DEDUP') ?: ''

process dedup_bam{
    conda "$DEDUPENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$DEDUPENV"
    cpus 4
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.endsWith("_dedup.bam"))              "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("_dedup.bam.bai") > 0)  "MAPPED/${COMBO}/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf("dedup.log") > 0)       "LOGS/${COMBO}/${CONDITION}/DEDUP/${file(filename).getName()}"
        else null
    }

    input:
    tuple val(sample_id), path(mapped_bam), path(ubam)
    path ref
        
    output:
    path "*_dedup.bam", emit: bam
    path "*_dedup.bam.bai", emit: bai
    path "*_dedup.log", emit: logs

    script:
    bams = mapped_bam
    outf = bams.getSimpleName()+"_dedup.bam"
    outl = bams.getSimpleName()+"_dedup.log"
    """
    mkdir -p TMP
    ref_dict=\$(basename $ref .gz).dict
    if [[ ! -f "\${ref_dict}" ]]; then samtools dict $ref -o \${ref_dict} >> $outl 2>&1; fi
    samtools sort -n -@ ${task.cpus} -o TMP/ubam_qn.bam $ubam >> $outl 2>&1
    samtools sort -n -@ ${task.cpus} -o TMP/mapped_qn.bam $bams >> $outl 2>&1
    samtools view -h TMP/mapped_qn.bam | awk 'BEGIN{FS=OFS="\t"} /^@/{print; next} {f=\$2+0; if (!and(f,256) && !and(f,2048)) {k=\$1":"(and(f,64)?1:0)":"(and(f,128)?1:0); if (seen[k]++) next} print}' | samtools view -b -o TMP/mapped_qn_primaryuniq.bam - >> $outl 2>&1
    $DEDUPBIN zipper --unmapped TMP/ubam_qn.bam --input TMP/mapped_qn_primaryuniq.bam --reference $ref --output TMP/zippered.bam --threads ${task.cpus} --compression-level 1 >> $outl 2>&1
    $DEDUPBIN sort --order template-coordinate --input TMP/zippered.bam --output TMP/sorted.bam --threads ${task.cpus} --compression-level 1 >> $outl 2>&1
    $DEDUPBIN dedup $DEDUPPARAMS --input TMP/sorted.bam --output TMP/dedup.bam >> $outl 2>&1
    $DEDUPBIN sort --order coordinate --input TMP/dedup.bam --output $outf --write-index --threads ${task.cpus} --compression-level 1 >> $outl 2>&1
    """
}

workflow DEDUPBAM{
    take:
    map
    mapi
    mapu
    mapui
    ubam

    main:
    mapped_ch = map.concat(mapu).map { b ->
        def n = file(b).getName()
        def key = n
            .replaceFirst(/_mapped_sorted_unique\.bam$/, '')
            .replaceFirst(/_mapped_sorted\.bam$/, '')
            .replaceFirst(/_R1_dedup_trimmed$/, '')
            .replaceFirst(/_dedup_trimmed$/, '')
            .replaceFirst(/_R1_trimmed$/, '')
            .replaceFirst(/_trimmed$/, '')
        tuple(key, b)
    }
    ubam_ch = ubam.map { u ->
        def n = file(u).getName()
        def key = n
            .replaceFirst(/_fgumi_extract\.bam$/, '')
            .replaceFirst(/_extracted\.bam$/, '')
            .replaceFirst(/_R1$/, '')
        tuple(key, u)
    }
    paired_ch = mapped_ch.combine(ubam_ch, by: 0).map { key, mb, ub -> tuple(key, mb, ub) }

    ref_ch = channel.value(file(REFERENCE))
    dedup_bam(paired_ch, ref_ch)

    emit:
    dedup = dedup_bam.out.bam
    dedupbai = dedup_bam.out.bai
    deduplog = dedup_bam.out.logs
}
