process make_rmd{
    conda "summary.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename == 'SUMMARY.html')         "REPORTS/SUMMARY/SUMMARY.html"                               
        else if (filename.indexOf("log") > 0)        "LOGS/REPORTS/SUMMARY/make_rmd.log"
    }

    input:
    path figs
    path tables

    output:
    path "*.html", emit: report
    path "log", emit: log

    script:
    """
    ln -f \"${projectDir}/../REPORTS/SUMMARY/summary.Rmd\" .;
    touch log;
    Rscript --vanilla -e \"rmarkdown::render('summary.Rmd', params=list(root='.'), output_file='SUMMARY.html', quiet=TRUE)\" 2> log
    """
}


workflow SUMMARY{
    take: collection

    main:

    png_ch =  Channel.fromPath("${projectDir}/../D{E,EU,AS,TU}/**/Figures/*.png")
    tab_ch =  Channel.fromPath("${projectDir}/../D{E,EU,AS,TU}/**/Tables/*.tsv.gz")
    //png_ch.subscribe {  println "PNG: $it"  }
    //tab_ch.subscribe {  println "TABLE: $it"  }

    make_rmd(png_ch.collect(), tab_ch.collect())
    
    emit:
    rmds = make_rmd.out.report
}

workflow{
    SUMMARY(dummy)
}