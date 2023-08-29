BINS = get_always('BINS')
DASENV = get_always('DASENV')
DASBIN = get_always('DASBIN')
DASREF = get_always('DASREF')
DASREFDIR = "${workflow.workDir}/../"+get_always('DASREFDIR')
DASANNO = get_always('DASANNO')
COUNTPARAMS = get_always('diego_DAS_params_COUNT') ?: ''
DASPARAMS = get_always('diego_DAS_params_DAS') ?: ''
DASSAMPLES = get_always('DASSAMPLES') ?: ''
DASGROUPS = get_always('DASREPS') ?: ''
DASREPS = get_always('DASREPS') ?: ''
DASCOMP = get_always('DASCOMP') ?: ''
DASCOMPS = get_always('DASCOMPS') ?: ''
PVAL = get_always('DASPVAL') ?: ''
LFC = get_always('DASLFC') ?: ''
PCOMBO = get_always('COMBO') ?: 'none'

COUNTBIN = 'featureCounts'
COUNTENV = 'countreads_de'

//DAS PROCESSES

process featurecount_diego{
    conda "$COUNTENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".counts.gz") > 0)     "DAS/${SCOMBO}/Featurecounts/${CONDITION}/${file(filename).getName()}"
        else if (filename.indexOf(".counts.summary") > 0)      "DAS/${SCOMBO}/Featurecounts/${CONDITION}/${file(filename).getName()}"             
        else if (filename.indexOf(".log") > 0)        "LOGS/DAS/${SCOMBO}/${file(filename).getSimpleName()}/featurecounts_diego_unique.log"
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

    output: 
    path "samplemap.txt", emit: smap
    path "groupings.txt", emit: groups
    path "log", emit: log

    script:
    """
    echo $DASSAMPLES|sed -e 's/:/\\t/' -e 's/\\|/\\n/' 1> samplemap.txt 2>> log && echo $DASGROUPS|sed -e 's/:/\\t/' -e 's/\\|/\\n/' 1> groupings.txt 2>> log
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
    python ${BINS}/Analysis/DAS/FeatureCounts2DIEGO.py $DASREPS --table junction_table.txt.gz --anno ANNOTATION.gz 2> log
    """
}

process create_contrast_files{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("contrast.txt") > 0 )      "DAS/${SCOMBO}/Tables/${COMBO}_contrast.txt"
        else if (filename == "log")      "LOGS/DAS/${SCOMBO}/${COMBO}_create_contrast_files.log"
    }

    input:
    path anno

    output: 
    path "*contrast.txt", emit: contrast
    path "log", emit: log

    script:
    """
    python ${BINS}/Analysis/DAS/diego_contrast_files.py -a <(zcat $anno) -b $DASCOMPS -c $DASCOMP -o . 2> log
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
    path tbl
    path contrast

    output:
    path "*.pdf", emit: dendrogram
    path "*.csv", emit: table
    path "log", emit: log

    script:    
    //outdir = "DAS"+File.separatorChar+"${SCOMBO}"
    bin = "${BINS}"+File.separatorChar+"${DASBIN}"
    outcomb = "DAS_DIEGO_"+"${COMBO}"
    """    
    arr=($contrast); for i in \${!arr[@]}; do basecond=\$(head -n 1 \${arr[\$i]} | awk \'{print \$1}\'); outcond=\$(echo \${basecond}|sed \'s/_[0-9]+//g\'); $bin -a <(zcat $tbl) -b \${arr[\$i]} -x \${basecond} -e -f ${outcomb}_\${outcond}_figure_dendrogram &> log; $bin -a <(zcat $tbl) -b \${arr[\$i]} -x \$basecond $DASPARAMS 1> ${outcomb}_\${outcond}_table_results.csv 2>> log;done
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
    set +o pipefail; arr=($tabs); for i in \${!arr[@]}; do a=\${arr[\$i]}; fn=\${a##*/}; if [[ -s \"\$a\" ]];then cat \$a|head -n1 > Sig_\$a; cat \$a| tail -n+2 |grep -v -w 'NA'|perl -F'\\t' -wlane 'next if (!\$F[10]);if (\$F[10] eq \"yes\") {print}' >> Sig_\$a &>> log; else touch \${orr[\$i]}; fi; done
    """
}

process convertPDF{
    conda "$DASENV"+".yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("dendrogram") > 0)      "DAS/${SCOMBO}/Figures/${file(filename).getName()}"                               
    }

    input:
    path pdf

    output:
    path "*.png", emit: png

    script:
    
    """
    for pdfile in $pdf ; do convert -verbose -density 500 -resize '800' \$pdfile \${pdfile%pdf}png; done
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
        else if (filename.indexOf("log") > 0)        "LOGS/DAS/filter_diego.log"
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

process collect_diego{
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

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())
    annofile = Channel.fromPath(DASANNO)
    
    featurecount_diego(annofile.combine(mapsamples_ch.collate(1)))
    create_samplemaps()
    prepare_junction_usage_matrix(create_samplemaps.out.smap, featurecount_diego.out.fc_cts.collect())
    create_contrast_files(prepare_junction_usage_matrix.output.anno)
    run_diego(prepare_junction_usage_matrix.out.jtab,create_contrast_files.out.contrast)
    filter_significant(run_diego.out.table)
    collect_diego(filter_significant.out.sigtbls.collect())
    convertPDF(run_diego.out.dendrogram)
    create_summary_snippet(run_diego.out.table.concat(filter_significant.out.sigtbls.concat(run_diego.out.dendrogram.concat(convertPDF.out.png))).collect())
    

    emit:
    tbls = run_diego.out.table
    sigtbls = filter_significant.out.sigtbls
    figs = run_diego.out.dendrogram
    snps = create_summary_snippet.out.snps
}