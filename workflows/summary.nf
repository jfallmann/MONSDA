process make_rmd{
    conda "summary.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/SUMMARY.html"                               
        else if (filename.indexOf("log") > 0)        "LOGS/REPORTS/SUMMARY/make_rmd.log"
    }

    input:
    path snippets
    path figs
    path tables

    output:
    path "*.html", emit: report
    path "log", emit: log

    script:
    inlist = snippets.toString()
    """
    touch log; Rscript --vanilla -e \"rmarkdown::render('$inlist', params=list(root='${workflow.workDir}/../'), output_file='SUMMARY.html', quiet=TRUE)\" 2> log
    """
}


workflow SUMMARY{ 
    take: collection

    main:

    sum_ch =  Channel.fromPath("${projectDir}/../REPORTS/SUMMARY/summary.Rmd")
    png_ch =  Channel.fromPath("${projectDir}/../{DE,DEU,DAS,DTU}/**/Figures/*.png")
    tab_ch =  Channel.fromPath("${projectDir}/../{DE,DEU,DAS,DTU}/**/Tables/*.tsv.gz")
    //png_ch.subscribe {  println "PNG: $it"  }
    //tab_ch.subscribe {  println "TABLE: $it"  }

    make_rmd(sum_ch, png_ch, tab_ch)
    
    emit:
    rmds = make_rmd.out.report
}

workflow{
    SUMMARY(dummy)
}