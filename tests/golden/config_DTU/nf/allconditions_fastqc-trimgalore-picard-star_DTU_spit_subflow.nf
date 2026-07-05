#!/usr/bin/env nextflow
nextflow.enable.dsl=2
def get_always(parameter){
    return params.containsKey(parameter) ? params[parameter] : null
}
params.gREFERENCE = "${workflow.workDir}/../"+get_always('REFERENCE')
params.gREFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
params.gBINS = { def BINS = get_always('BINS'); BINS = get_always('BINS'); return BINS }()
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
params.gDTUENV = get_always('DTUENV')
params.gDTUBIN = get_always('DTUBIN')
params.gDTUREF = get_always('DTUREF')
params.gDTUREFDIR = "${workflow.workDir}/../"+get_always('DTUREFDIR')
params.gDTUANNO = get_always('DTUANNO')
params.gDTUDECOY = get_always('DTUDECOY')
params.gDTUIDX = get_always('DTUIDX')
params.gDTUUIDX = get_always('DTUUIDX')
params.gDTUUIDXNAME = get_always('DTUUIDXNAME')+'.idx'
params.gIDXPARAMS = get_always('spit_DTU_params_INDEX') ?: ''
params.gCOUNTPARAMS = get_always('spit_DTU_params_COUNT') ?: ''
params.gDTUPARAMS = get_always('spit_DTU_params_DTU') ?: ''
params.gRUNTERMINUS = get_always('spit_DTU_params_TERMINUS') ?: ''
params.gDTUREPS = get_always('DTUREPS') ?: ''
params.gDTUCOMP = get_always('DTUCOMP') ?: ''
params.gDTUCOMPS = get_always('DTUCOMPS') ?: ''
params.gPVAL = get_always('DTUPVAL') ?: ''
params.gLFC = get_always('DTULFC') ?: ''
params.gPCOMBO = get_always('params.gCOMBO') ?: 'none'
params.gCOUNTBIN = 'salmon'
params.gCOUNTENV = 'salmon'
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

process salmon_quant{
    conda "<REPO>/envs/${params.gCOUNTENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gCOUNTENV}"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${params.gSCOMBO}/salmon/${params.gCONDITION}/DTU/${file(filename).getName()}"
        else                                    "DTU/${params.gSCOMBO}/salmon/${params.gCONDITION}/"+"${filename.replaceAll(/trimmed./,"")}"
    }
    input:
    path reads
    output:
    path "*.gz", emit: counts
    path "*.log", emit: logs
    script:
    idx = reads[0]
    if (params.gPAIRED == 'paired'){
        if (params.gSTRANDED == 'fr' || params.gSTRANDED == 'ISF'){
            stranded = '-l ISF'
        }else if (params.gSTRANDED == 'rf' || params.gSTRANDED == 'ISR'){
            stranded = '-l ISR'
        }else{
            stranded = '-l IU'
        }
        rs = reads[1..2].sort { a,b -> a[0] <=> b[0] == 0 ? (a[1..-1] as int) <=> (b[1..-1] as int) : a[0] <=> b[0] }
        r1 = rs[0]
        r2 = rs[1]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        lf = "salmon_"+fn+".log"
        of = fn+"/quant.sf"
        oz = fn+"/quant.sf.gz"
        ol = fn+"_counts.gz"
        """
        ${params.gCOUNTBIN} ${params.gCOUNTPARAMS} quant -p ${task.cpus} -i $idx $stranded -o $fn -1 $r1 -2 $r2 &>> $lf && gzip $of && ln -fs $oz $ol
        """
    }else{
        if (params.gSTRANDED == 'fr' || params.gSTRANDED == 'SF'){
            stranded = '-l SF'
        }else if (params.gSTRANDED == 'rf' || params.gSTRANDED == 'SR'){
            stranded = '-l SR'
        }else{
            stranded = '-l U'
        }
        read = reads[1]
        fn = file(reads[1]).getSimpleName().replaceAll(/\Q_trimmed\E/,"")
        lf = "salmon_"+fn+".log"
        of = fn+"/quant.sf"
        oz = fn+"/quant.sf.gz"
        ol = fn+"_counts.gz"
        """
        ${params.gCOUNTBIN} ${params.gCOUNTPARAMS} quant -p ${task.cpus} -i $idx $stranded -o $fn -r $read &>> $lf && gzip $of && ln -fs $oz $ol
        """
    }
}
process create_summary_snippet{
    conda "<REPO>/envs/${params.gDTUENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gDTUENV}"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${params.gSCOMBO}.Rmd"                               
        else if (filename.indexOf("log") > 0)        "LOGS/DTU/create_summary_snippet.log"
    }
    input:
    path de
    output:
    path "*.Rmd", emit: snps
    path "log", emit: log
    script:
    inlist = de.toString()
    """
    touch log; python3 ${params.gBINS}/Analysis/RmdCreator.py --files $inlist --output out.Rmd --env ${params.gDTUENV} --loglevel DEBUG 2>> log
    """
}
process terminus_collapse{
    conda "<REPO>/envs/terminus.yaml"
    container "oras://docker.io/jfallmann/monsda:terminus"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${params.gSCOMBO}/terminus/${params.gCONDITION}/DTU/${file(filename).getName()}"
        else                                    "DTU/${params.gSCOMBO}/terminus/${file(filename).getName()}"
    }
    input:
    path counts
    output:
    path "*.gz", emit: counts
    path "*.log", emit: logs
    script:
    lf = "terminus_collapse.log"
    """
    mkdir -p termdirs; for f in *_counts.gz; do b=\${f%_counts.gz}; mkdir -p termdirs/\$b; zcat \$f > termdirs/\$b/quant.sf; done 2> $lf; dirs=\$(ls -d termdirs/*); for d in \$dirs; do terminus group ${params.gRUNTERMINUS} -d \$d -o termout &>> $lf; done; terminus collapse -d \$dirs -o termout &>> $lf; for d in \$dirs; do b=\$(basename \$d); gzip -c termout/\$b/quant.sf > \${b}_counts.gz; done 2>> $lf
    """
}
process salmon_idx{
    conda "<REPO>/envs/${params.gCOUNTENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gCOUNTENV}"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)    "LOGS/${params.gCOMBO}/${params.gCONDITION}/DTU/spit_index.log"
        else if (filename == "spit.idx")           "${params.gDTUIDX}"
        else                                          "${params.gDTUUIDX}"
    }
    input:
    path genome
    output:
     path "${params.gDTUUIDXNAME}", emit: idx
    path "*.log", emit: idxlog
    path "*.idx", emit: tmpidx
    script:    
    gen =  genome.getName()
    if (params.gDTUDECOY){
        decoy = "-d "+"${params.gDTUDECOY}" 
    }else{
        decoy = ''
    }
    """
    ${params.gCOUNTBIN} index ${params.gIDXPARAMS} $decoy -p ${task.cpus} -t $gen -i ${params.gDTUUIDXNAME} &> index.log && ln -fs ${params.gDTUUIDXNAME} spit.idx
    """
}
process prepare_dtu_annotation{
    conda "<REPO>/envs/${params.gDTUENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gDTUENV}"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".gz") > 0)       "DTU/${params.gSCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/DTU/${params.gSCOMBO}/featurecount_spit_annotation.log"
    }
    output:
    path "*.gz", emit: anno
    path "*.log", emit: log
    script:     
    ca = params.gCOMBO+"_ANNOTATION.gz"
    ol = "create_DTU_table.log"
    """
    mkdir -p TMP; ${params.gBINS}/Analysis/build_DTU_table.py ${params.gDTUREPS} --anno $ca --loglevel DEBUG --nextflow 2>> $ol
    """
}
process run_spit{
    conda "<REPO>/envs/${params.gDTUENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gDTUENV}"
    cpus params.gTHREADS
	cache 'lenient'
    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DTU/${params.gSCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("_figure") > 0)      "DTU/${params.gSCOMBO}/Figures/${file(filename).getName()}" 
        else if (filename.indexOf("SESSION") > 0)      "DTU/${params.gSCOMBO}/${file(filename).getName()}"                     
        else if (filename.indexOf("log") > 0)        "LOGS/DTU/${params.gSCOMBO}/run_spit.log"
    }
    input:
    path counts
    path anno
    path ref
    output:
    path "*_table*", emit: tbls
    path "*_figure*", emit: figs, optional:true
    path "*SESSION.gz", emit: session
    path "log", emit: log
    script:    
    outdir = "DTU"+File.separatorChar+"${params.gSCOMBO}"
    bin = "${params.gBINS}"+File.separatorChar+"${params.gDTUBIN}"
    comp = "${params.gDTUCOMP}".split(':')[0]
    dparams = "'${params.gDTUPARAMS}'"
    """
    mkdir -p Figures Tables
    python3 $bin $anno $ref . ${params.gDTUCOMP} ${params.gPCOMBO} ${task.cpus} $dparams &> log ; ln -f Tables/* . && touch Figures/dummy() && ln -f Figures/* .
    """
}
process collect_spit{
    conda "<REPO>/envs/${params.gDTUENV}"+".yaml"
    container "oras://docker.io/jfallmann/monsda:"+"${params.gDTUENV}"
    cpus params.gTHREADS
	cache 'lenient'
    input:
    path de
    script:
    """
    echo "$de DONE"
    """
}
workflow DTU{ 
    take: collection
    main:
    TRIMSAMPLES = params.gLONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${params.gCOMBO}/"+element+"_{R2,R1}*.fastq.gz"
        }
    trimsamples_ch =  Channel.fromPath(TRIMSAMPLES.sort())
    annofile = Channel.fromPath(params.gDTUANNO)
    checkidx = file(params.gDTUIDX)
    if (checkidx.exists()){
        idxfile = Channel.fromPath(params.gDTUUIDX)
        if (params.gPAIRED == 'paired'){
            salmon_quant(idxfile.combine(trimsamples_ch.collate(2)))
        } else{
            salmon_quant(idxfile.combine(trimsamples_ch.collate(1)))
        }        
    }
    else{
        genomefile = Channel.fromPath(params.gDTUREF)
        salmon_idx(genomefile)
        if (params.gPAIRED == 'paired'){
            salmon_quant(salmon_idx.out.idx.combine(trimsamples_ch.collate(2)))
        } else{
            salmon_quant(salmon_idx.out.idx.combine(trimsamples_ch.collate(1)))
        }
    }
    prepare_dtu_annotation()
    if (params.gRUNTERMINUS?.length() > 0){
        terminus_collapse(salmon_quant.out.counts.collect())
        dtucounts = terminus_collapse.out.counts
    }
    else{
        dtucounts = salmon_quant.out.counts
    }
    run_spit(dtucounts.collect(), prepare_dtu_annotation.out.anno, annofile)
    create_summary_snippet(run_spit.out.tbls.concat(run_spit.out.figs.concat(run_spit.out.session)).collect())
    collect_spit(run_spit.out.tbls.collect().concat(create_summary_snippet.out.snps.collect()))
    emit:
    tbls = run_spit.out.tbls
    figs = run_spit.out.figs
    snps = create_summary_snippet.out.snps
}




workflow {
    DTU(dummy())
}

