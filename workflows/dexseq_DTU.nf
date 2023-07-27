BINS = get_always('BINS')
DTUENV = get_always('DTUENV')
DTUBIN = get_always('DTUBIN')
DTUREF = get_always('DTUREF')
DTUREFDIR = "${workflow.workDir}/../"+get_always('DTUREFDIR')
DTUANNO = get_always('DTUANNO')
DTUDECOY = get_always('DTUDECOY')
DTUIDX = get_always('DTUIDX')
DTUUIDX = get_always('DTUUIDX')
DTUUIDXNAME = get_always('DTUUIDXNAME')+'.idx'
IDXPARAMS = get_always('dexseq_DTU_params_INDEX') ?: ''
COUNTPARAMS = get_always('dexseq_DTU_params_COUNT') ?: ''
DTUPARAMS = get_always('dexseq_DTU_params_DTU') ?: ''
DTUREPS = get_always('DTUREPS') ?: ''
DTUCOMP = get_always('DTUCOMP') ?: ''
DTUCOMPS = get_always('DTUCOMPS') ?: ''
PVAL = get_always('DTUPVAL') ?: ''
LFC = get_always('DTULFC') ?: ''
PCOMBO = get_always('COMBO') ?: 'none'

COUNTBIN = 'salmon'
COUNTENV = 'salmon'

//DTU PROCESSES

process salmon_idx{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/DTU/dexseq_index.log"
        else if (filename == "dexseq.idx")            "$DTUIDX"
        else                                          "$DTUUIDX"
    }

    input:
    path genome

    output:
     path "$DTUUIDXNAME", emit: idx
    path "*.log", emit: idxlog
    path "*.idx", emit: tmpidx

    script:    
    gen =  genome.getName()
    if (DTUDECOY){
        decoy = "-d "+"${DTUDECOY}" 
    }else{
        decoy = ''
    }
    """
    $COUNTBIN index $IDXPARAMS $decoy -p ${task.cpus} -t $gen -i $DTUUIDXNAME &> index.log && ln -fs $DTUUIDXNAME dexseq.idx
    """

}

process salmon_quant{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${SCOMBO}/salmon/${CONDITION}/DTU/${file(filename).getName()}"
        else                                    "DTU/${SCOMBO}/salmon/${CONDITION}/"+"${filename.replaceAll(/trimmed./,"")}"
    }

    input:
    path reads

    output:
    path "*.gz", emit: counts
    path "*.log", emit: logs

    script:

    idx = reads[0]
    if (PAIRED == 'paired'){
        if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '-l ISF'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
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
        $COUNTBIN $COUNTPARAMS quant -p ${task.cpus} -i $idx $stranded -o $fn -1 $r1 -2 $r2 &>> $lf && gzip $of && ln -fs $oz $ol
        """
    }else{
        if (STRANDED == 'fr' || STRANDED == 'SF'){
            stranded = '-l SF'
        }else if (STRANDED == 'rf' || STRANDED == 'SR'){
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
        $COUNTBIN $COUNTPARAMS quant -p ${task.cpus} -i $idx $stranded -o $fn -r $read &>> $lf && gzip $of && ln -fs $oz $ol
        """
    }
}

process prepare_dtu_annotation{
    conda "$DTUENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".gz") > 0)       "DTU/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf(".log") > 0)        "LOGS/DTU/${SCOMBO}/featurecount_dexseq_annotation.log"
    }

    output:
    path "*.gz", emit: anno
    path "*.log", emit: log

    script:     
    ca = COMBO+"_ANNOTATION.gz"
    ol = "create_DTU_table.log"
    """
    mkdir -p TMP; $BINS/Analysis/build_DTU_table.py $DTUREPS --anno $ca --loglevel DEBUG --nextflow 2>> $ol
    """
}

process run_dexseq{
    conda "$DTUENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DTU/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("_figure") > 0)      "DTU/${SCOMBO}/Figures/${file(filename).getName()}" 
        else if (filename.indexOf(".html") > 0)      "DTU/${SCOMBO}/DEXSeqReport_${COMBO}_${DTUCOMP}/${file(filename).getName()}"
        else if (filename.indexOf("SESSION") > 0)      "DTU/${SCOMBO}/${file(filename).getName()}"                     
        else if (filename.indexOf("log") > 0)        "LOGS/DTU/${SCOMBO}/run_dexseq.log"
    }

    input:path counts
    path anno
    path ref

    output:
    path "*_table*", emit: tbls
    path "*_figure*", emit: figs, optional:true
    path "*.html", emit: html, optional:true
    path "*SESSION.gz", emit: session
    path "log", emit: log

    script:    
    outdir = "DTU"+File.separatorChar+"${SCOMBO}"
    bin = "${BINS}"+File.separatorChar+"${DTUBIN}"
    comp = "${DTUCOMP}".split(':')[0]
    dparams = "'${DTUPARAMS}'"

    """
    mkdir -p Figures Tables
    Rscript --no-environ --no-restore --no-save $bin $anno $ref . $DTUCOMP $PCOMBO ${task.cpus} $dparams &> log ; ln -f Tables/* . && touch Figures/dummy && ln -f Figures/* .
    """
}

process create_summary_snippet{
    conda "$DTUENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${SCOMBO}.Rmd"                               
        else if (filename.indexOf("log") > 0)        "LOGS/DTU/create_summary_snippet.log"
    }

    input:
    path de

    output:
    path "*.Rmd", emit: snps
    path "log", emit: log

    script:
    inlist = de.toString()
    // inlist = de.toList()  // { $workflow.projectDir += "$it.code,"  }
    """
    touch log; python3 $BINS/Analysis/RmdCreator.py --files $inlist --output out.Rmd --env $DTUENV --loglevel DEBUG 2>> log
    """
}

process collect_dexseq{
    conda "$DTUENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

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
    
    TRIMSAMPLES = LONGSAMPLES.collect{
            element -> return "${workflow.workDir}/../TRIMMED_FASTQ/${COMBO}/"+element+"_{R2,R1}*.fastq.gz"
        }

    trimsamples_ch =  Channel.fromPath(TRIMSAMPLES.sort())
    annofile = Channel.fromPath(DTUANNO)
    checkidx = file(DTUIDX)
    
    if (checkidx.exists()){
        idxfile = Channel.fromPath(DTUUIDX)
        if (PAIRED == 'paired'){
            salmon_quant(idxfile.combine(trimsamples_ch.collate(2)))
        } else{
            salmon_quant(idxfile.combine(trimsamples_ch.collate(1)))
        }        
    }
    else{
        genomefile = Channel.fromPath(DTUREF)
        salmon_idx(genomefile)
        if (PAIRED == 'paired'){
            salmon_quant(salmon_idx.out.idx.combine(trimsamples_ch.collate(2)))
        } else{
            salmon_quant(salmon_idx.out.idx.combine(trimsamples_ch.collate(1)))
        }
    }

    prepare_dtu_annotation()
    run_dexseq(salmon_quant.out.counts.collect(), prepare_dtu_annotation.out.anno, annofile)
    create_summary_snippet(run_dexseq.out.tbls.concat(run_dexseq.out.figs.concat(run_dexseq.out.session)).collect())
    collect_dexseq(run_dexseq.out.tbls.collect().concat(create_summary_snippet.out.snps.collect()))

    emit:
    tbls = run_dexseq.out.tbls
    figs = run_dexseq.out.figs
    snps = create_summary_snippet.out.snps
}