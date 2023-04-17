BINS = get_always('BINS')
DTUENV = get_always('DTUENV')
DTUBIN = get_always('DTUBIN')
DTUREF = get_always('DTUREF')
DTUREFDIR = "${workflow.workDir}/../"+get_always('DTUREFDIR')
DTUANNO = get_always('DTUANNO')
DTUIDX = get_always('DTUIDX')
DTUUIDX = get_always('DTUUIDX')
DTUUIDXNAME = get_always('DTUUIDXNAME')
IDXPARAMS = get_always('drimseq_DTU_params_INDEX') ?: ''
COUNTPARAMS = get_always('drimseq_DTU_params_COUNT') ?: ''
DTUPARAMS = get_always('drimseq_DTU_params_DTU') ?: ''
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
        if (filename == "salmon.idx")            "$COUNTUIDX"
        else if (filename.indexOf(".log") >0)    "LOGS/${COMBO}/${CONDITION}/COUNTING/salmon_index.log"
    }

    input:
    path genome

    output:
    path "salmon.idx", emit: idx

    script:    
    gen =  genome.getName()
    """
    $COUNTBIN index $IDXPARAMS -p $THREADS -t $gen -i $DTUUIDXNAME &> index.log && ln -fs $DTUUIDXNAME salmon.idx
    """

}

process salmon_quant{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".sf.gz") >0)            "COUNTS/${SCOMBO}/${CONDITION}/"+"${filename.replaceAll(/trimmed./,"")}"
        else if (filename.indexOf(".log") >0)               "LOGS/${SCOMBO}/${CONDITION}/COUNTING/${file(filename).getName()}"
        else null
    }

    input:
    path reads

    output:
    path "*.sf.gz", emit: counts
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
        r1 = reads[1]
        r2 = reads[2]
        fn = file(r1).getSimpleName().replaceAll(/\Q_R1_trimmed\E/,"")
        lf = "salmon_"+fn+".log"
        of = fn+"/quant.sf"
        oz = fn+"/quant.sf.gz"
        ol = fn+"_counts.sf.gz"
        """
        $COUNTBIN $COUNTPARAMS quant -p $THREADS -i $idx $stranded -o $fn -1 $r1 -2 $r2 &>> $lf && gzip $of && mv -f $oz $ol
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
        ol = fn+"_counts.sf.gz"
        """
        $COUNTBIN $COUNTPARAMS quant -p $THREADS -i $idx $stranded -o $fn -r $read &>> $lf && gzip $of && mv -f $oz $ol
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
        else if (filename.indexOf(".log") > 0)        "LOGS/DTU/${SCOMBO}/featurecount_drimseq_annotation.log"
    }

    output:
    path "*.gz", emit: anno
    path "*.log", emit: log

    script:     
    ca = COMBO+"_ANNOTATION.gz"
    ol = "create_DTU_table.log"
    """
    mkdir -p TMP; $BINS/Analysis/build_DTU_table.py $DTUREPS --anno $ca --loglevel DEBUG 2>> $ol
    """
}

process run_drimseq{
    conda "$DTUENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DTU/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("_figure") > 0)      "DTU/${SCOMBO}/Figures/${file(filename).getName()}" 
        else if (filename.indexOf(".html") > 0)      "DTU/${SCOMBO}/drimseqReport_${COMBO}_${DTUCOMP}/${file(filename).getName()}"
        else if (filename.indexOf("SESSION") > 0)      "DTU/${SCOMBO}/${file(filename).getName()}"                     
        else if (filename.indexOf("log") > 0)        "LOGS/DTU/${SCOMBO}/run_drimseq.log"
    }

    input:
    path anno
    path ref

    output:
    path "*_table*", emit: tbls
    path "*_figure*", emit: figs
    path "*.html", emit: html
    path "*SESSION.gz", emit: session
    path "log", emit: log

    script:    
    outdir = "DTU"+File.separatorChar+"${SCOMBO}"
    bin = "${BINS}"+File.separatorChar+"${DTUBIN}"

    """
    mkdir -p Figures Tables drimseqReport_${COMBO}_${DTUCOMP}
    Rscript --no-environ --no-restore --no-save $bin $anno $ref . $PCOMBO $DTUCOMP $THREADS $DTUPARAMS 2> log && mv Tables/* . && mv Figures/* . && mv drimseqReport_*/* .
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

process collect_drimseq{
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

    trimsamples_ch =  Channel.fromPath(TRIMSAMPLES)
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
    run_drimseq(prepare_dtu_annotation.out.anno, annofile)
    create_summary_snippet(run_drimseq.out.tbls.concat(run_drimseq.out.figs.concat(run_drimseq.out.session)).collect())
    collect_drimseq(run_drimseq.out.tbls.collect().concat(create_summary_snippet.out.snps.collect()))

    emit:
    tbls = run_drimseq.out.tbls
    figs = run_drimseq.out.figs
    snps = create_summary_snippet.out.snps
}