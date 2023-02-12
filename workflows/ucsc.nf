BINS = get_always('BINS')
TRACKSENV = get_always('TRACKSENV')
TRACKSBIN = get_always('TRACKSBIN')
REF = get_always('REFERENCE')
REFDIR = get_always('REFDIR')
ANNO = get_always('ANNO')
TRACKSPARAMS = get_always('ucsc_TRACKS_params_UCSC') ?: ''

TRACKBIN = 'ucsc'
TRACKENV = 'ucsc'

include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"

process BamToBed{
    conda "bedtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bam
    
    output:
    path "*.bed.gz", emit: bed
    path "*.log", emit: log
    
    script: 
    fn = file(bam).getSimpleName()
    fo = fn+'.bed.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    if (STRANDED == 'rf' || STRANDED == 'ISR'){
        """
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/1\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=$THREADS -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }else{
        """
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=$THREADS -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }
}


process BedToBedg{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_ucscbedtobedgraph.log"
    }

    input:
    path bed
    path fai
    path sizes

    output:
    path "*fw.bedg.gz", emit: bedgf
    path "*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    fn = file(bed).getSimpleName()
    fw = fn+'.fw.bedg.gz'
    fr = fn+'.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | sort --parallel=$THREADS -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fw 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes |sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}


process NormalizeBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_ucscnormalizebedgraph.log"
    }

    input:
    path bedgf
    path bedgr

    output:
    path "*.fw.norm.bedg.gz", emit: bedgf
    path "*.re.norm.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.norm.fw.bedg.gz'
    fr = fn+'.norm.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgf|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\t\",@F[0..\$#F-1]),\"\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgr|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\t\",@F[0..\$#F-1]),\"\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S $sortmem -T TMP -t\$'\t' -k1,1 -k2,2n|gzip > $bedgr 2> $ol; else gzip < /dev/null > $bedgr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process BedgToTRACKS{
    conda "ucsc.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/${SCOMBO}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_bedgtoucsc.log"
    }

    input:
    path bedgf
    path bedgr
    path sizes

    output:
    path "*.fw.bw", emit: bwf
    path "*.re.bw", emit: bwr
    path "*.log", emit: log

    script: 
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.bw'
    fr = fn+'.re.bw'
    ol = fn+".log"
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\n\' __)\" ]] ;then zcat $bedgf > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\n\' __)\" ]] ;then zcat $bedgr > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process GenerateTrack{
    conda "base.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".txt") > 0)      "TRACKS/${file(filename).getName()}"
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${file(filename).getName()}_track.log"
    }

    input:
    path bwf
    path bwr

    output:
    path "*.txt", emit: trackdb
    path "*.log", emit: log

    script: 
    fn = file(bwf).getSimpleName()
    ol = fn+".log"
    src='TRACKS'+os.sep+$SETS.replace(os.sep, '_')
    
    """
    echo -e \"$bwf\n$bwr\"|python3 $BINS/Analysis/GenerateTrackDb.py -i $SETS -e 1 -f STDIN -u '' -g $REFDIR $TRACKSPARAMS 2> $ol
    """
}

workflow TRACKS{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted_*.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES)
    mapsamples_ch.subscribe {  println "MAP: $it \t COMBO: ${COMBO} SCOMBO: ${SCOMBO} LONG: ${LONGSAMPLES}"  }
    genomefile = Channel.fromPath(REF)

    UnzipGenome(genomefile)
    BamToBed(mapsamples_ch.collate(1))
    BedToBedg(BamToBed.out.bed, UnzipGenome.out.index, UnzipGenome.out.chromsize)
    NormalizeBedg(BedToBedg.out.bedgf, BedToBedg.out.bedgr)
    BedgToTRACKS(NormalizeBedg.out.bedgf, NormalizeBedg.out.bedgr, UnzipGenome.out.chromsize)
    GenerateTrack(BedgToTRACKS.out.bwf.collect(), BedgToTRACKS.out.bwr.collect())

    emit:
    trackdb = GenerateTrack.out.trackdb
}
