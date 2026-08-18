FUSENV = get_always('FUSIONSENV')
FUSBIN = get_always('FUSIONSBIN')
FUSREF = get_always('FUSIONSREF')
FUSREFDIR = "${workflow.workDir}/../"+get_always('FUSIONSREFDIR')
FUSANNO = get_always('FUSIONSANNO')
FUSLIB = get_always('FUSIONSLIB')

FUSPARAMS = get_always('starfusion_params_FUSION') ?: ''
FUSBUILD = get_always('starfusion_params_BUILD') ?: ''

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
        else      "FUSIONS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path fls

    output:
    path "${fn}_starfusion/*", emit: fusions
    path "log", emit: log

    script:
    ref = fls[0]
    anno = fls[1]
    junction = fls[2]        
    fn = file(junction).getSimpleName()
    od = fn+"_starfusion"
    ol = fn+".log"
    
    """
    if [[ -f \"$FUSLIB/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIB\"; elif [[ -f \"$FUSLIB/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIB/ctat_genome_lib_build_dir\"; else CTAT=CTAT; ( mkdir -p \$CTAT && zcat ${ref} > ctat_ref.fa && zcat ${anno} > ctat_ref.gtf && prep_genome_lib.pl --genome_fa ctat_ref.fa --gtf ctat_ref.gtf --output_dir \$CTAT --CPU ${task.cpus} $FUSBUILD ) &> log; fi; if [[ -s \"${junction}\" ]]; then ctat_chr=\$(grep -v '^#' \$CTAT/ref_annot.gtf 2>/dev/null | head -1 | cut -f1 | grep -c '^chr' || true); junc_chr=\$(awk '!/^#/ && \$1!=\"chr_donorA\"{print \$1; exit}' ${junction} | grep -c '^chr' || true); if [[ \"\$ctat_chr\" == \"1\" && \"\$junc_chr\" == \"0\" ]]; then echo \"Adding chr prefix to junction to match CTAT lib\" >> log; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {if(\$1!~/^chr/)\$1=\"chr\"\$1; if(\$4!~/^chr/)\$4=\"chr\"\$4; if(\$1==\"chrMT\")\$1=\"chrM\"; if(\$4==\"chrMT\")\$4=\"chrM\"; print}' ${junction} > ${fn}.norm.junction; elif [[ \"\$ctat_chr\" == \"0\" && \"\$junc_chr\" == \"1\" ]]; then echo \"Stripping chr prefix from junction to match CTAT lib\" >> log; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {sub(/^chr/,\"\",\$1); sub(/^chr/,\"\",\$4); if(\$1==\"M\")\$1=\"MT\"; if(\$4==\"M\")\$4=\"MT\"; print}' ${junction} > ${fn}.norm.junction; else cp ${junction} ${fn}.norm.junction; fi; $FUSBIN --genome_lib_dir \$CTAT -J ${fn}.norm.junction --output_dir ${od} --CPU ${task.cpus} $FUSPARAMS &>> log; else mkdir -p ${od}; echo \"File ${junction} empty, no chimeric STAR output found\" >> log; fi; mkdir -p ${od}; touch ${od}/star-fusion.fusion_predictions.tsv ${od}/star-fusion.fusion_predictions.abridged.tsv
    """
}

workflow FUSIONS{ 
    take: collection

    main:

    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"*.Chimeric.out.junction"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())  
    annofile = Channel.fromPath(FUSANNO)
    genomefile = Channel.fromPath(FUSREF)

    starfusion(genomefile.combine(annofile.combine(mapsamples_ch.collate(1))))

    emit:
    fusions = starfusion.out.fusions
    logs = starfusion.out.log
}
