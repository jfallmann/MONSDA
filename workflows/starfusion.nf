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
        if (filename.indexOf("_fusions") > 0)      "FUSIONS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}"        
        else if (filename.indexOf(".log") > 0)        "LOGS/${SCOMBO}/${CONDITION}/${file(filename).getSimpleName()}"
    }

    input:
    path fls

    output:
    path "*_fusions", emit: fusions
    path "log", emit: log

    script:
    ref = fls[0]
    anno = fls[1]
    junction = fls[2]        
    fn = file(junction).getSimpleName()
    of = fn+"_fusions"
    ol = fn+".log"
    
    """
    if [[ -f \"$FUSLIB/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIB\"; elif [[ -f \"$FUSLIB/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then CTAT=\"$FUSLIB/ctat_genome_lib_build_dir\"; else CTAT=CTAT; ( mkdir -p \$CTAT && zcat ${ref} > ctat_ref.fa && zcat ${anno} > ctat_ref.gtf && prep_genome_lib.pl --genome_fa ctat_ref.fa --gtf ctat_ref.gtf --output_dir \$CTAT --CPU ${task.cpus} $FUSBUILD ) &> log; fi; if [[ -s \"${junction}\" ]]; then $FUSBIN --genome_lib_dir \$CTAT -J ${junction} --output_dir TMP --CPU ${task.cpus} $FUSPARAMS &>> log && cp TMP/star-fusion.fusion_predictions.tsv ${of} 2>> log; else echo \"File ${junction} empty, no chimeric STAR output found\" >> log; fi; touch ${of}
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
