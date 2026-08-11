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
params.gFUSENV = get_always('FUSIONSENV')
params.gFUSBIN = get_always('FUSIONSBIN')
params.gFUSREF = get_always('FUSIONSREF')
params.gFUSREFDIR = "${workflow.workDir}/../"+get_always('FUSIONSREFDIR')
params.gFUSANNO = get_always('FUSIONSANNO')
params.gFUSLIB = get_always('FUSIONSLIB')
params.gFUSPARAMS = get_always('starfusion_params_FUSION') ?: ''
params.gFUSBUILD = get_always('starfusion_params_BUILD') ?: ''
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

process starfusion{
    conda "<REPO>/envs/${params.gFUSENV}"+".yaml"
    container "oras://ghcr.io/jfallmann/monsda:"+"${params.gFUSENV}"+"-VERSION"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".log") > 0)        "LOGS/${params.gSCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
        else      "FUSIONS/${params.gSCOMBO}/${params.gCONDITION}/${file(filename).getName()}"
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
    if [[ -f \"${params.gFUSLIB}/ref_annot.gtf\" ]]; then CTAT=\"${params.gFUSLIB}\"; elif [[ -f \"${params.gFUSLIB}/ctat_genome_lib_build_dir/ref_annot.gtf\" ]]; then CTAT=\"${params.gFUSLIB}/ctat_genome_lib_build_dir\"; else CTAT=CTAT; ( mkdir -p \$CTAT && zcat ${ref} > ctat_ref.fa && zcat ${anno} > ctat_ref.gtf && prep_genome_lib.pl --genome_fa ctat_ref.fa --gtf ctat_ref.gtf --output_dir \$CTAT --CPU ${task.cpus} ${params.gFUSBUILD} ) &> log; fi; if [[ -s \"${junction}\" ]]; then ctat_chr=\$(grep -v '^#' \$CTAT/ref_annot.gtf 2>/dev/null | head -1 | cut -f1 | grep -c '^chr' || true); junc_chr=\$(awk '!/^#/ && \$1!=\"chr_donorA\"{print \$1; exit}' ${junction} | grep -c '^chr' || true); if [[ \"\$ctat_chr\" == \"1\" && \"\$junc_chr\" == \"0\" ]]; then echo \"Adding chr prefix to junction to match CTAT lib\" >> log; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {if(\$1!~/^chr/)\$1=\"chr\"\$1; if(\$4!~/^chr/)\$4=\"chr\"\$4; if(\$1==\"chrMT\")\$1=\"chrM\"; if(\$4==\"chrMT\")\$4=\"chrM\"; print}' ${junction} > ${fn}.norm.junction; elif [[ \"\$ctat_chr\" == \"0\" && \"\$junc_chr\" == \"1\" ]]; then echo \"Stripping chr prefix from junction to match CTAT lib\" >> log; awk 'BEGIN{OFS=\"\\t\"} /^#/{print;next} \$1==\"chr_donorA\"{print;next} {sub(/^chr/,\"\",\$1); sub(/^chr/,\"\",\$4); if(\$1==\"M\")\$1=\"MT\"; if(\$4==\"M\")\$4=\"MT\"; print}' ${junction} > ${fn}.norm.junction; else cp ${junction} ${fn}.norm.junction; fi; ${params.gFUSBIN} --genome_lib_dir \$CTAT -J ${fn}.norm.junction --output_dir ${od} --CPU ${task.cpus} ${params.gFUSPARAMS} &>> log; else mkdir -p ${od}; echo \"File ${junction} empty, no chimeric STAR output found\" >> log; fi; mkdir -p ${od}; touch ${od}/star-fusion.fusion_predictions.tsv ${od}/star-fusion.fusion_predictions.abridged.tsv
    """
}
workflow FUSIONS{ 
    take: collection
    main:
    MAPPEDSAMPLES = params.gLONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${params.gCOMBO}/"+element+"*.Chimeric.out.junction"
    }
    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())  
    annofile = Channel.fromPath(params.gFUSANNO)
    genomefile = Channel.fromPath(params.gFUSREF)
    starfusion(genomefile.combine(annofile.combine(mapsamples_ch.collate(1))))
    emit:
    fusions = starfusion.out.fusions
    logs = starfusion.out.log
}




workflow {
    FUSIONS(dummy())
}

