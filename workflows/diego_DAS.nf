BINS = get_always('BINS')
DASENV = get_always('DASENV')
DASBIN = get_always('DASBIN')
DASREF = get_always('DASREF')
DASREFDIR = get_always('DASREFDIR')
DASANNO = get_always('DASANNO')
COUNTPARAMS = get_always('edger_DAS_params_COUNT') ?: ''
DASPARAMS = get_always('edger_DAS_params_DAS') ?: ''
DASREPS = get_always('DASREPS') ?: ''
DASCOMP = get_always('DASCOMP') ?: ''
DASCOMPS = get_always('DASCOMPS') ?: ''
PVAL = get_always('DASPVAL') ?: ''
LFC = get_always('DASLFC') ?: ''
PCOMBO = get_always('COMBO') ?: 'none'

COUNTBIN = 'featureCounts'
COUNTENV = 'countreads_de'

//DAS PROCESSES

process featurecount_edger{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".count") > 0)      "DAS/${SCOMBO}/Featurecounts/${file(filename).getSimpleName()}.counts.gz"                
        else if (filename.indexOf(".log") > 0)        "LOGS/DAS/${SCOMBO}/${file(filename).getSimpleName()}/featurecounts_edger_unique.log"
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
    mkdir -p TMP; $COUNTBIN -T $THREADS $COUNTPARAMS $pair $stranded -a <(zcat $anno) -o tmpcts $reads 2> $ol && head -n2 tmpcts |gzip > $oc && export LC_ALL=C; tail -n+3 tmpcts|sort --parallel=$THREADS -S $sortmem -T TMP -k1,1 -k2,2n -k3,3n -u |gzip >> $oc 2>> $ol && mv tmpcts.summary $os
    """
}


process create_samplemaps{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "samplemap.txt")      "DAS/${SCOMBO}/Tables/${COMBO}_samplemap.txt"
        else if (filename == "groupings.txt")      "DAS/${SCOMBO}/Tables/${COMBO}_grouping.txt"
        else if (filename == "log")      "LOGS/DAS/${SCOMBO}/${COMBO}_create_samplemaps.log"
    }

    input:
    //path '*.count*'// from reads
    path samples
    path groups

    output: 
    path "samplemap.txt", emit: smap
    path "groupings.txt", emit: groups
    path "log", emit: log

    script:
    """
    echo $samples 1> samplemap.txt 2>> log && echo $groups 1> groupings.txt 2>> log
    """
}


process prepare_junction_usage_matrix{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "junction_table.txt.gz")      "DAS/${SCOMBO}/Tables/${COMBO}_junction_table_dexdas.txt.gz"
        else if (filename == "ANNOTATION.gz")      "DAS/${SCOMBO}/Tables/${COMBO}_ANNOTATION.gz"
        else if (filename == "log")      "LOGS/DAS/${SCOMBO}/${COMBO}_junction_usage_matrix.log"
    }

    input:
    //path '*.count*'// from reads
    path smap
    path cts

    output: 
    path "junction_table.txt.gz", emit: jtab
    path "*ANNOTATION.gz", emit: anno
    path "log", emit: log

    script:
    """
    ${BINS}/Analysis/DAS/FeatureCounts2DIEGO.py $DASREPS --table junction_table.txt.gz --anno ANNOTATION.gz 2> log
    """
}

process create_contrast_files{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == "contrast.txt")      "DAS/${SCOMBO}/Tables/${COMBO}_contrast.txt"
        else if (filename == "log")      "LOGS/DAS/${SCOMBO}/${COMBO}_create_contrast_files.log"
    }

    input:
    //path '*.count*'// from reads
    path jtab

    output: 
    path "contrast.txt", emit: contrast
    path "*ANNOTATION.gz", emit: anno
    path "log", emit: log

    script:
    """
    ${BINS}/Analysis/DAS/diego_contrast_files.py  -a <(zcat $jtab) -b $DASCOMP -c $DASCOMPS -o . 2> log
    """
}

process run_diego{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("csv") > 0)      "DAS/${SCOMBO}/Tables/${file(filename).getName()}"                
        else if (filename.indexOf("dendrogram") > 0)      "DAS/${SCOMBO}/Figures/${file(filename).getName()}"                
        else if (filename.indexOf("log") > 0)        "LOGS/DAS/${SCOMBO}_${COMBO}_${COMPSTR}/run_diego.log"
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
    outdir = "DAS"+File.separatorChar+"${SCOMBO}"
    bin = "${BINS}"+File.separatorChar+"${DASBIN}"
    """
    mkdir -p Figures Tables
    Rscript --no-environ --no-restore --no-save $bin $anno $cts $deanno . $DASCOMP $PCOMBO $THREADS $DASPARAMS 2> log && mv Tables/* . && mv Figures/* .
    """
}

process filter_significant{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_table") > 0)      "DAS/${SCOMBO}/Tables/${file(filename).getName()}"                                
        else if (filename.indexOf("log") > 0)        "LOGS/DAS/filter_deseq2.log"
    }

    input:
    path tabs

    output:
    path "*_table*", emit: sigtbls
    path "log", emit: log

    script:  
    """
    set +o pipefail; for i in $tabs; do if [[ -s \"\${i}\" ]];then zcat \${i}| head -n1 |gzip > Sig_\${i};cp -f Sig_\${i} SigUP_\${i}; cp -f Sig_\${i} SigDOWN_\${i}; zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] <= -$LFC ||\$F[3] >= $LFC) ){{print}}' |gzip >> Sig_\${i} && zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] >= $LFC) ){{print}}' |gzip >> SigUP_\${i} && zcat \${i}| tail -n+2 |grep -v -w 'NA'|perl -F\'\\t\' -wlane 'next if (!\$F[6] || !\$F[3]);if (\$F[6] < $PVAL && (\$F[3] <= -$LFC) ){{print}}' |gzip >> SigDOWN_\${i}; else touch Sig_\${i} SigUP\${i} SigDOWN_\${i}; fi;done 2> log
    """
}

process create_summary_snippet{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${SCOMBO}.Rmd"                               
        else if (filename.indexOf("log") > 0)        "LOGS/DAS/filter_edger.log"
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
    touch log; python3 $BINS/Analysis/RmdCreator.py --files $inlist --output out.Rmd --env $DASENV --loglevel DEBUG 2>> log
    """
}

process collect_edger{
    conda "$DASENV"+".yaml"
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

workflow DAS{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted_unique.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES)
    mapsamples_ch.subscribe {  println "MAP: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }
    annofile = Channel.fromPath(DASANNO)
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