BINS = get_always('BINS')
DEENV = get_always('DEENV')
DEBIN = get_always('DEBIN')
DEREF = get_always('DEREF')
DEREFDIR = "${workflow.workDir}/../"+get_always('DEREFDIR')
DEANNO = get_always('DEANNO')
COUNTPARAMS = get_always('edger_DE_params_COUNT') ?: ''
DEPARAMS = get_always('edger_DE_params_DE') ?: ''
DEREPS = get_always('DEREPS') ?: ''
DECOMP = get_always('DECOMP') ?: ''
DECOMPS = get_always('DECOMPS') ?: ''
PVAL = get_always('DEPVAL') ?: ''
LFC = get_always('DELFC') ?: ''
PCOMBO = get_always('COMBO') ?: 'none'

COUNTBIN = 'featureCounts'
COUNTENV = 'countreads_de'

//DE PROCESSES

process featurecount_edger{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".counts.gz") > 0)      "DE/${SCOMBO}/Featurecounts/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".counts.summary") > 0)      "DE/${SCOMBO}/Featurecounts/${CONDITION}/${file(filename).getName()}"               
        else if (filename.indexOf(".log") > 0)        "LOGS/DE/${SCOMBO}/${file(filename).getSimpleName()}/featurecounts_edger_unique.log"
    }

    input:
    path fls

    output:
    path "*.counts.gz", emit: fc_cts
    path "*.summary", emit: fc_summary
    path "*.log", emit: fc_log

    script: 
    anno = fls[0]
    reads = fls[1]       
    fn = file(reads).getSimpleName()
    oc = fn+".counts.gz"
    os = fn+".counts.summary"
    ol = fn+".log"
    sortmem = '30%'
    if (PAIRED == 'paired'){
        pair = "-p"
    }
    else{
        pair= ""
    }
    if (STRANDED == 'fr' || STRANDED == 'ISF'){
            stranded = '-s 1'
        }else if (STRANDED == 'rf' || STRANDED == 'ISR'){
            stranded = '-s 2'
        }else{
            stranded = ''
    }
    """
    mkdir -p TMP; $COUNTBIN -T ${task.cpus} $COUNTPARAMS $pair $stranded -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=${task.cpus} -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc 2>> $ol && mv tmpcts.summary $os
    """
}

process prepare_count_table{
    conda "$DEENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "COUNTS.gz")      "DE/${SCOMBO}/Tables/${COMBO}_COUNTS.gz"
        else if (filename == "ANNOTATION.gz")      "DE/${SCOMBO}/Tables/${COMBO}_ANNOTATION.gz"
        else if (filename == "SampleDict.gz")      "DE/${SCOMBO}/Tables/${COMBO}_SampleDict.gz"
        else if (filename == "log")      "LOGS/DE/${SCOMBO}/${COMBO}_prepare_count_table.log"
    }

    input:
    //path '*.count*'// from reads
    path reps

    output: 
    path "*COUNTS.gz", emit: counts
    path "*ANNOTATION.gz", emit: anno
    path "*SampleDict.gz", emit: sdict
    path "log", emit: log

    script:
    """
    ${BINS}/Analysis/build_count_table.py $DEREPS --table COUNTS.gz --anno ANNOTATION.gz --nextflow 2> log
    """
}

process run_edger{
    conda "$DEENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DE/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("_figure") > 0)      "DE/${SCOMBO}/Figures/${file(filename).getName()}"                
        else if (filename.indexOf("SESSION") > 0)      "DE/${SCOMBO}/${file(filename).getName()}"                     
        else if (filename.indexOf("log") > 0)        "LOGS/DE/${SCOMBO}/run_edger.log"
    }

    input:
    //path '*.count*'// from reads
    path cts
    path anno
    path deanno

    output:
    path "*_table*", emit: tbls
    path "*_figure*", emit: figs
    path "*SESSION.gz", emit: session
    path "log", emit: log

    script:    
    outdir = "DE"+File.separatorChar+"${SCOMBO}"
    bin = "${BINS}"+File.separatorChar+"${DEBIN}"
    """
    mkdir -p Figures Tables
    Rscript --no-environ --no-restore --no-save $bin $anno $cts $deanno . $DECOMP $PCOMBO ${task.cpus} $DEPARAMS 2> log && mv Tables/* . && mv Figures/* .
    """
}

process filter_significant{
    conda "$DEENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DE/${SCOMBO}/Tables/${file(filename).getName()}"                                
        else if (filename.indexOf("log") > 0)        "LOGS/DE/filter_deseq2.log"
    }

    input:
    path tabs

    output:
    path "*_table*", emit: sigtbls
    path "log", emit: log

    script:  
    """
    set +o pipefail; for i in $tabs; do if [[ -s \"\${i}\" ]];then zcat \${i}| head -n1 |gzip > Sig_\${i};cp -f Sig_\${i} SigUP_\${i}; cp -f Sig_\${i} SigDOWN_\${i}; zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F'\\t' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] <= -$LFC ||\$F[3] >= $LFC) ){{print}}' |gzip >> Sig_\${i} && zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F'\\t' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] >= $LFC) ){{print}}' |gzip >> SigUP_\${i} && zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F'\\t' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] <= -$LFC) ){{print}}' |gzip >> SigDOWN_\${i}; else touch Sig_\${i} SigUP\${i} SigDOWN_\${i}; fi;done 2> log
    """
}

process create_summary_snippet{
    conda "$DEENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${SCOMBO}.Rmd"                               
        else if (filename.indexOf("log") > 0)        "LOGS/DE/filter_edger.log"
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
    touch log; python3 $BINS/Analysis/RmdCreator.py --files $inlist --output out.Rmd --env $DEENV --loglevel DEBUG 2>> log
    """
}

process collect_edger{
    conda "$DEENV"+".yaml"
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

workflow DE{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted_unique.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())
    //mapsamples_ch.subscribe {  println "MAP: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }
    annofile = Channel.fromPath(DEANNO)
    //annofile.subscribe {  println "ANNO: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }

    featurecount_edger(annofile.combine(mapsamples_ch.collate(1)))
    prepare_count_table(featurecount_edger.out.fc_cts.collect())
    run_edger(prepare_count_table.out.counts, prepare_count_table.out.anno, annofile)
    filter_significant(run_edger.out.tbls)
    create_summary_snippet(run_edger.out.tbls.concat(run_edger.out.figs.concat(run_edger.out.session)).collect())
    collect_edger(filter_significant.out.sigtbls.collect())

    emit:
    tbls = run_edger.out.tbls
    sigtbls = filter_significant.out.sigtbls
    figs = run_edger.out.figs
    snps = create_summary_snippet.out.snps
}