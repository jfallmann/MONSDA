//DTU shared processes

process salmon_quant{
    conda "$COUNTENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$COUNTENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${SCOMBO}/${CONDITION}/DTU/salmon/${file(filename).getName()}"
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

process create_summary_snippet{
    conda "$DTUENV"+".yaml"
    container "oras://jfallmann/monsda:"+"$DTUENV"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".Rmd") > 0)         "REPORTS/SUMMARY/RmdSnippets/${SCOMBO}.Rmd"                               
        else if (filename.indexOf("log") > 0)        "LOGS/${SCOMBO}/DTU/salmon/create_summary_snippet.log"
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

process terminus_collapse{
    conda "terminus.yaml"
    container "oras://jfallmann/monsda:terminus"
    cpus THREADS
	cache 'lenient'

    publishDir "${workflow.workDir}/../" , mode: 'copyNoFollow',
    saveAs: {filename ->
        if (filename.indexOf(".log") >0)        "LOGS/${SCOMBO}/${CONDITION}/DTU/terminus/${file(filename).getName()}"
        else                                    "DTU/${SCOMBO}/terminus/${file(filename).getName()}"
    }

    input:
    path counts

    output:
    path "*.gz", emit: counts
    path "*.log", emit: logs

    script:
    lf = "terminus_collapse.log"
    """
    mkdir -p termdirs; for f in *_counts.gz; do b=\${f%_counts.gz}; mkdir -p termdirs/\$b; zcat \$f > termdirs/\$b/quant.sf; done 2> $lf; dirs=\$(ls -d termdirs/*); for d in \$dirs; do terminus group $RUNTERMINUS -d \$d -o termout &>> $lf; done; terminus collapse -d \$dirs -o termout &>> $lf; for d in \$dirs; do b=\$(basename \$d); gzip -c termout/\$b/quant.sf > \${b}_counts.gz; done 2>> $lf
    """
}
