BINS = get_always('BINS')
PEAKSENV = get_always('PEAKSENV')
PEAKSBIN = get_always('PEAKSBIN')
REF = get_always('REFERENCE')
REFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
SETS = get_always('SETS')
IP = get_always('IP')
SOFTPARAMS = get_always('peaks_PEAKS_params_SOFTCLIP') ?: ''
PREPARAMS = get_always('peaks_PEAKS_params_PREPROCESS') ?: ''
PEAKSPARAMS = get_always('peaks_PEAKS_params_FINDPEAKS') ?: ''

include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"

process BamToBed{
    conda "bedtools.yaml"
    cpus THREADS
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bam2bed.log"
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
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/1\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }else{
        """
        bedtools bamtobed -split -i $bam | sed 's/ /_/g' | perl -wl -a -F'\\t' -n -e '\$F[0] =~ s/\\s/_/g;if(\$F[3]=~/\\/2\$/){{if (\$F[5] eq \"+\"){{\$F[5] = \"-\"}}elsif(\$F[5] eq \"-\"){{\$F[5] = \"+\"}}}} print join(\"\\t\",@F[0..\$#F])' | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fo 2> $ol
        """
    }
}

process ExtendBed{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "BED/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bam2bed.log"
    }

    input:
    path bedf

    output:
    path "*d.bed.gz", emit: bedext
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    sizes = bedf[1]

    fn = file(bed).getSimpleName()
    ol = fn+".log"
    sortmem = '30%'

    if (IP == 'iCLIP'){
        of = fn+'_extended.bed.gz'    
        opt = '-u 1'
        
    } else if (IP == 'revCLIP'){
        of = fn+'_revtrimmed.bed.gz'    
        opt = '-d 1'
    }

    """        
    $BINS/Universal/ExtendBed.pl $opt -b $bed -o $of -g $sizes 2> $ol
    """    
}

process BedToBedg{
    conda "bedtools.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "TRACKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_bedtobedgraph.log"
    }

    input:
    path bedf

    output:
    path "*.bedg.gz", emit: bedg
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    fai = bedf[1]
    sizes = bedf[2]

    fn = file(bed).getSimpleName()
    of = fn+'.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'

    """
    export LC_ALL=C; export LC_COLLATE=C; bedtools genomecov -i $bed -bg -split -strand + -g $sizes | perl -wlane 'print join(\"\\t\",@F[0..2],\".\",\$F[3],\"+\")' > tosrt 2> $ol && bedtools genomecov -i $bed -bg -split -strand - -g $sizes | perl -wlane 'print join(\"\\t\",@F[0..2],\".\",\$F[3],\"-\")' >> tosrt 2>> $ol && cat tosrt | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $of 2>> $ol
    """
}

process PreprocessPeaks{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/prepeak_${file(filename).getName()}.log"
    }

    input:
    path bed

    output:
    path "*_prepeak*.bed.gz", emit: prepeak
    path "*.log", emit: log

    script: 
    of = file(bed).getSimpleName().replaceAll(/\Q_mapped\E/,"_prepeak")+".bed.gz"
    ol = file(bed).getSimpleName()+".log"
    sortmem = '30%'

    """        
    perl $BINS/Analysis/PreprocessPeaks.pl -p <(zcat $bed) $PREPARAMS | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k3,3n -k2,2n -k6,6 |gzip > $of 2> $ol
    """    
}

process FindPeaks{
    conda "$PEAKSENV"+".yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bed.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/findpeaks_${file(filename).getName()}.log"
    }

    input:
    path bed

    output:
    path "*_peak*.bed.gz", emit: peak
    path "*.log", emit: log

    script: 
    of = file(bed).getName().replaceAll(/\Q_prepeak\E/,"_peak")
    ol = file(bed).getSimpleName()+".log"
    sortmem = '30%'

    """  
    export LC_ALL=C; if [[ -n \"\$(zcat $bed | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then $PEAKSBIN $PEAKSPARAMS <(zcat $bed|sort -t\$\'\\t\' -k1,1 -k3,3n -k2,2n -k6,6) 2> $ol|tail -n+2| sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$\'\\t\' -k1,1 -k2,2n |grep -v \'nan\'| gzip > $of 2>> $ol; else gzip < /dev/null > $of; echo \"File $bed empty\" >> $ol; fi
    """    
}


process PeakToBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_peaktobedgraph.log"
    }

    input:
    path bedf

    output:
    path "*fw.bedg.gz", emit: bedgf
    path "*re.bedg.gz", emit: bedgr
    path "*.log", emit: log

    script: 
    bed = bedf[0]
    sizes = bedf[1]

    fn = file(bed).getSimpleName()
    fw = fn+'.fw.bedg.gz'
    fr = fn+'.re.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'

    """
    perl $BINS/Universal/Bed2Bedgraph.pl -f <(zcat $bed) -c $sizes -p peak -x tmp.fw.gz -y tmp.re.gz -a track 2>> $ol && zcat tmp.fw.gz | sort --parallel=${task.cpus} -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2>> $ol && zcat tmp.re.gz |sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fr 2>> $ol
    """
}


process NormalizeBedg{
    conda "perl.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".norm.bedg.gz") > 0)      "PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_ucscpeaknormalizebedgraph.log"
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
    fw = fn+'.fw.norm.bedg.gz'
    fr = fn+'.re.norm.bedg.gz'
    ol = fn+".log"
    sortmem = '30%'
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgf|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgf)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n |gzip > $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then scale=\$(bc <<< \"scale=6;\$(zcat $bedgr|cut -f4|perl -wne '{\$x+=\$_;}END{if (\$x == 0){\$x=1} print \$x}')/1000000\") perl -wlane '\$sc=\$ENV{scale};print join(\"\\t\",@F[0..\$#F-1]),\"\\t\",\$F[-1]/\$sc' <(zcat $bedgr)| sort -S $sortmem -T TMP -t\$'\\t' -k1,1 -k2,2n|gzip > $fr 2> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process PeakToTRACKS{
    conda "ucsc.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".bw") > 0)      "TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"                
        else if (filename == ".log")        "LOGS/TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}_peaktoucsc.log"
    }

    input:
    path bedgf
    path rest
    
    output:
    path "*.fw.bw", emit: bwf
    path "*.re.bw", emit: bwr
    path "*.log", emit: log

    script: 
    bedgr = rest[0]
    sizes = rest[1]
    fn = file(bedgf).getSimpleName()
    fw = fn+'.fw.bw'
    fr = fn+'.re.bw'
    ol = fn+".log"
    
    """
    export LC_ALL=C; if [[ -n \"\$(zcat $bedgf | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgf > tmp && bedGraphToBigWig tmp $sizes $fw 2> $ol; else gzip < /dev/null > $fw; echo \"File $bedgf empty\" >> $ol; fi && if [[ -n \"\$(zcat $bedgr | head -c 1 | tr \'\\0\\n\' __)\" ]] ;then zcat $bedgr > tmp && bedGraphToBigWig tmp $sizes $fr 2>> $ol; else gzip < /dev/null > $fr; echo \"File $bedgr empty\" >> $ol; fi
    """
}

process GenerateTrack{
    conda "base.yaml"
    cpus 1
	cache 'lenient'
    //validExitStatus 0,1

    publishDir "${workflow.workDir}/../" , mode: 'link',
    saveAs: {filename ->
        if (filename.indexOf(".txt") > 0)      "TRACKS/PEAKS/${file(filename).getName()}"
        else if (filename == ".log")        "LOGS/TRACKS/PEAKS/${SCOMBO}/${CONDITION}/${file(filename).getName()}"
    }

    input:
    path bwf
    path bwr

    output:
    path "*.txt", emit: trackdb
    path "*.log", emit: log

    script: 
    uid= SETS.replace(File.separator, "_")
    ol = uid+"_GenerateTrack_peaks.log"
    opt = '-n Peaks_'+"$PEAKSENV"+' -s peaks -l TRACKS_peaks_'+"$PEAKSENV"+' -b TRACKS_'+"$PEAKSENV"
    """
    mkdir -p LOGS;touch LOGS/MONSDA.log; bf=($bwf); br=($bwr); blen=\${#bf[@]}; for i in \"\${!bf[@]}\";do f=\${bf[\$i]}; r=\${br[\$i]}; echo -e \"\$f\\n\$r\"|python3 $BINS/Analysis/GenerateTrackDb.py -i $uid -e 1 -f STDIN -u \"TRACKS/$SETS\" -g $REFDIR $opt 2>> $ol;done
    """
}


workflow PEAKS{ 
    take: collection

    main:
    
    MAPPEDSAMPLES = LONGSAMPLES.collect{
        element -> return "${workflow.workDir}/../MAPPED/${COMBO}/"+element+"_mapped_sorted_*.bam"
    }

    mapsamples_ch = Channel.fromPath(MAPPEDSAMPLES.sort())
    genomefile = Channel.fromPath(REF)

    UnzipGenome(genomefile)
    BamToBed(mapsamples_ch.collate(1))    
    if (IP == 'iCLIP' || IP == 'revCLIP'){
        ExtendBed(BamToBed.out.bed.combine(UnzipGenome.out.chromsize))
        BedToBedg(ExtendBed.out.bedext.combine(UnzipGenome.out.index.combine(UnzipGenome.out.chromsize)))
    } else {
        BedToBedg(BamToBed.out.bed.combine(UnzipGenome.out.index.combine(UnzipGenome.out.chromsize)))    
    }
    PreprocessPeaks(BedToBedg.out.bedg)
    FindPeaks(PreprocessPeaks.out.prepeak)
    PeakToBedg(FindPeaks.out.peak.combine(UnzipGenome.out.chromsize))
    NormalizeBedg(PeakToBedg.out.bedgf.collate(1), PeakToBedg.out.bedgr.collate(1))
    PeakToTRACKS(NormalizeBedg.out.bedgf.collate(1), NormalizeBedg.out.bedgr.collate(1).combine(UnzipGenome.out.chromsize))
    GenerateTrack(PeakToTRACKS.out.bwf.collect(), PeakToTRACKS.out.bwr.collect())

    emit:
    trackdb = GenerateTrack.out.trackdb
}
