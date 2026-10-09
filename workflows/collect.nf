process collect_stuff{
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->        
        "LOGS/${COMBO}/${CONDITION}/COLLECT/collect/${file(filename).getName()}"
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
