process collect_stuff{
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->        
        "LOGS/COLLECT/$COMBO$CONDITION/${file(filename).getName()}"
    }
    input:
    path check


    output:
    path "collect.txt", emit: done

    script:
    """
    echo "$check successful!" > collect.txt
    """
}

workflow COLLECT{
    take:
    whatever

    main:
    
    collect_stuff(whatever.collect())

    
}
