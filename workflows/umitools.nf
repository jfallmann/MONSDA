DEDUPENV=get_always('DEDUPENV')
DEDUPBIN=get_always('DEDUPBIN')

WHITELISTPARAMS = get_always('umitools_params_WHITELIST') ?: ''
EXTRACTPARAMS = get_always('umitools_params_EXTRACT') ?: ''

process whitelist{
    conda "$DEDUPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_whitelist") > 0)         "DEDUP_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}_whitelist"
        else if (filename.indexOf("log") > 0)           "LOGS/$COMBO$CONDITION/DEDUP/dedup_whitelist.log"
        else null
    }

    input:
    path samples
        
    output:
    path "*_whitelist", emit: wl

    script:    
    if (WHITELISTPARAMS == ''){    
        outf = samples[0].getSimpleName().replace("_R1","")+"_dummy_whitelist"
        """
        touch $outf
        """
    } else {
        if (PAIRED == 'paired'){
            r1 = samples[0]
            r2 = samples[1]
            outf = samples[0].getSimpleName().replace("_R1","")+"_whitelist"
            """
                mkdir tmp && $DEDUPBIN whitelist $WHITELISTPARAMS --temp-dir tmp --log=wl.log --stdin=$r1 --read2-in=$r2 --stdout=$outf
            """
        }
        else{
            outf = samples.getSimpleName()+"_whitelist"
            """
                mkdir tmp && $DEDUPBIN whitelist $WHITELISTPARAMS --temp-dir tmp --log=wl.log --stdin=$samples --stdout=$outf
            """
        }
    }
}

process extract_fq{
    conda "$DEDUPENV"+".yaml"
    cpus THREADS
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf("_dedup.fastq.gz") > 0)      "DEDUP_FASTQ/$COMBO$CONDITION/${file(filename).getSimpleName()}.fastq.gz"
        else if (filename.indexOf("log") > 0)             "LOGS/$COMBO$CONDITION/DEDUP/dedup_extract.log"
        else null
    }

    input:
    path wl
    path samples
        
    output:
    path "*_dedup.fastq.gz", emit: extract
    path "ex.log", emit: logs

    script:
    if (PAIRED == 'paired'){
        r1 = samples[0]
        r2 = samples[1]
        outf = samples[0].getSimpleName()+"_dedup.fastq.gz"
        outf2 = samples[1].getSimpleName()+"_dedup.fastq.gz"
        if (wl ==~ /\*dummy_whitelist/){
            """
                mkdir tmp && $DEDUPBIN extract $EXTRACTPARAMS --temp-dir tmp --log=ex.log --stdin=$r1 --read2-in=$r2 --stdout=$outf --read2-out=$outf2
            """
        }
        else{
            """
                mkdir tmp && $DEDUPBIN extract $EXTRACTPARAMS --whitelist=$wl --temp-dir tmp --log=ex.log --stdin=$r1 --read2-in=$r2 --stdout=$outf --read2-out=$outf2
            """
        }
    }
    else{
        outf = samples.getSimpleName()+"_dedup.fastq.gz"
        if (wl ==~ 'dummy_whitelist'){
            """
                mkdir tmp && $DEDUPBIN extract $EXTRACTPARAMS --temp-dir tmp --log=ex.log --stdin=$samples --stdout=$outf
            """
        }
        else{        
            """
                mkdir tmp && $DEDUPBIN extract $EXTRACTPARAMS --whitelist=$wl --temp-dir tmp --log=ex.log --stdin=$samples --stdout=$outf
            """
        }
    }
}

workflow DEDUPEXTRACT{
    take: 
    collection

    main:
    //SAMPLE CHANNELS
    if ( PREDEDUP == 'enabled' ){ 
        if (PAIRED == 'paired'){                
            whitelist(samples_ch.collate( 2 ))
            extract_fq(whitelist.out.wl, samples_ch.collate( 2 ))
        } else{                
            whitelist(samples_ch.collate( 1 ))
            extract_fq(whitelist.out.wl, samples_ch.collate( 1 ))
        }
    }else{
        if (PAIRED == 'paired'){
            whitelist(collection.collate(2))
            extract_fq(whitelist.out.wl, collection.collate( 2 ))
        } else{
            whitelist(collection.collate( 1 ))
            extract_fq(whitelist.out.wl, collection.collate( 1 ))
        }
    }

    emit:    
    extract = extract_fq.out.extract
    logs = extract_fq.out.logs
}