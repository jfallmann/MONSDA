FUSENV = get_always('FUSIONSENV')
FUSBIN = get_always('FUSIONSBIN')
FUSREF = get_always('FUSIONSREF')
FUSREFDIR = "${workflow.workDir}/../"+get_always('FUSIONSREFDIR')
FUSANNO = get_always('FUSIONSANNO')
FUSLIB = get_always('FUSIONSLIB')
FUSLIBDIR = FUSLIB.startsWith('/') ? FUSLIB : "${workflow.workDir}/../"+FUSLIB

FUSPARAMS = get_always('starfusion_params_FUSION') ?: ''
FUSBUILD = get_always('starfusion_params_BUILD') ?: ''
FUSFASTQ = (get_always('starfusion_params_FASTQ') ?: '').toString().toLowerCase() in ['1', 'true', 'yes']

//FUSIONS PROCESSES

process starfusion{
    conda "$FUSENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$FUSENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/FUSIONS/starfusion/${file(filename).getName()}"
        else      "FUSIONS/${SCOMBO}/${CONDITION}/"+filename.toString().replace('_starfusion/', '/')
    }

    input:
    path fls

    output:
    path "${fn}_starfusion/*", emit: fusions
    path "*.log", emit: log

    script:
    ref = fls[0]
    anno = fls[1]
    junction = fls[2]        
    fn = file(junction).getSimpleName()
    od = fn+"_starfusion"
    ol = fn+".log"
    
    """
    if [[ -f \"$FUSLIBDIR/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIBDIR\"; elif [[ -f \"$FUSLIBDIR/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIBDIR/ctat_genome_lib_build_dir\"; else CTAT=CTAT; ( mkdir -p \$CTAT && zcat ${ref} > ctat_ref.fa && zcat ${anno} > ctat_ref.gtf && prep_genome_lib.pl --genome_fa ctat_ref.fa --gtf ctat_ref.gtf --output_dir \$CTAT --CPU ${task.cpus} $FUSBUILD ) &> ${ol}; fi; if [[ -s \"${junction}\" ]]; then ctat_chr=\$(grep -v '^#' \$CTAT/ref_annot.gtf 2>/dev/null | head -1 | cut -f1 | grep -c '^chr' || true); junc_chr=\$(awk '!/^#/ && \$1!=\"chr_donorA\"{print \$1; exit}' ${junction} | grep -c '^chr' || true); if [[ \"\$ctat_chr\" == \"1\" && \"\$junc_chr\" == \"0\" ]]; then echo \"Adding chr prefix to junction to match CTAT lib\" >> ${ol}; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {if(\$1!~/^chr/)\$1=\"chr\"\$1; if(\$4!~/^chr/)\$4=\"chr\"\$4; if(\$1==\"chrMT\")\$1=\"chrM\"; if(\$4==\"chrMT\")\$4=\"chrM\"; print}' ${junction} > ${fn}.norm.junction; elif [[ \"\$ctat_chr\" == \"0\" && \"\$junc_chr\" == \"1\" ]]; then echo \"Stripping chr prefix from junction to match CTAT lib\" >> ${ol}; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {sub(/^chr/,\"\",\$1); sub(/^chr/,\"\",\$4); if(\$1==\"M\")\$1=\"MT\"; if(\$4==\"M\")\$4=\"MT\"; print}' ${junction} > ${fn}.norm.junction; else cp ${junction} ${fn}.norm.junction; fi; $FUSBIN --genome_lib_dir \$CTAT -J ${fn}.norm.junction --output_dir ${od} --CPU ${task.cpus} $FUSPARAMS &>> ${ol}; else mkdir -p ${od}; echo \"File ${junction} empty, no chimeric STAR output found\" >> ${ol}; fi; mkdir -p ${od}; touch ${od}/star-fusion.fusion_predictions.tsv ${od}/star-fusion.fusion_predictions.abridged.tsv
    """
}

process starfusion_fastq{
    conda "$FUSENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$FUSENV"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/FUSIONS/starfusion/${file(filename).getName()}"
        else      "FUSIONS/${SCOMBO}/${CONDITION}/"+filename.toString().replace('_starfusion/', '/')
    }

    input:
    path fls

    output:
    path "${fn}_starfusion/*", emit: fusions
    path "*.log", emit: log

    script:
    ref = fls[0]
    anno = fls[1]
    if (PAIRED == 'paired'){
        rs = fls[2..3].sort{ it.getName() }
        reads = "--left_fq "+rs[0]+" --right_fq "+rs[1]
        fn = file(rs[0]).getSimpleName().replaceAll(/_R[12]_trimmed/, "")
    } else{
        reads = "--left_fq "+fls[2]
        fn = file(fls[2]).getSimpleName().replaceAll(/_trimmed/, "")
    }
    od = fn+"_starfusion"
    ol = fn+".log"

    """
    if [[ -f \"$FUSLIBDIR/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIBDIR\"; elif [[ -f \"$FUSLIBDIR/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIBDIR/ctat_genome_lib_build_dir\"; else CTAT=CTAT; ( mkdir -p \$CTAT && zcat ${ref} > ctat_ref.fa && zcat ${anno} > ctat_ref.gtf && prep_genome_lib.pl --genome_fa ctat_ref.fa --gtf ctat_ref.gtf --output_dir \$CTAT --CPU ${task.cpus} $FUSBUILD ) &> ${ol}; fi; $FUSBIN --genome_lib_dir \$CTAT ${reads} --output_dir ${od} --CPU ${task.cpus} $FUSPARAMS &>> ${ol}; mkdir -p ${od}; touch ${od}/star-fusion.fusion_predictions.tsv ${od}/star-fusion.fusion_predictions.abridged.tsv
    """
}

workflow FUSIONS{ 
    take: collection

    main:

    annofile = Channel.fromPath(FUSANNO)
    genomefile = Channel.fromPath(FUSREF)

    if (FUSFASTQ){
        if (PAIRED == 'paired'){
            TRIMSAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"_{R2,R1}_trimmed.fastq.gz"
            }
        } else{
            TRIMSAMPLES = LONGSAMPLES.collect{
                element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"_trimmed.fastq.gz"
            }
        }

        trimsamples_ch = Channel.fromPath(TRIMSAMPLES.sort())

        if (PAIRED == 'paired'){
            starfusion_fastq(genomefile.combine(annofile.combine(trimsamples_ch.collate(2))))
        } else{
            starfusion_fastq(genomefile.combine(annofile.combine(trimsamples_ch.collate(1))))
        }
        outfusions = starfusion_fastq.out.fusions
        outlogs = starfusion_fastq.out.log
    }
    else{
        MAPPEDSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"*.Chimeric.out.junction"
        }

        mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())

        starfusion(genomefile.combine(annofile.combine(mapsamples_ch.collate(1))))
        outfusions = starfusion.out.fusions
        outlogs = starfusion.out.log
    }

    emit:
    fusions = outfusions
    logs = outlogs
}
