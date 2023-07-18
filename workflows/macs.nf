BINS = get_always('BINS')
PEAKSENV = get_always('PEAKSENV')
PEAKSBIN = get_always('PEAKSBIN')
REF = get_always('REFERENCE')
REFDIR = "${workflow.workDir}/../"+get_always('REFDIR')
SETS = get_always('SETS')
IP = get_always('IP')
PREPARAMS = get_always('macs_params_PREPROCESS') ?: ''
PEAKSPARAMS = get_always('macs_params_FINDPEAKS') ?: ''
PEAKSAMPLES = get_always('PEAKSAMPLES')

include { UnzipGenome; UnzipGenome_no_us } from "manipulate_genome.nf"

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
    path bam

    output:
    path "*_peak.bed.gz", emit: peak
    path "*.log", emit: log

    script: 
    bf = bam[0]
    bc = bam[1]
    fn = file(bf).getSimpleName()
    of = fn+"_peak.bed"
    oz = fn+"_peak.bed.gz"
    ol = fn+".log"
    sortmem = '30%'
    if (PAIRED == 'paired' && bam.indexOf("unique") == 0){
        mapmode = 'BAMPE'
    }else{
        mapmode = 'BAM'
    }
    
    """      
    set +o pipefail; export LC_ALL=C; if [[ -n \"\$(samtools view $bf | head -c 1 | tr '\\0\\n' __)\" ]] ;then $PEAKBIN callpeak -t $bf -c $bc --outdir . -n $of -f $mapmode $PEAKSPARAMS 2> $ol && gzip *_peaks.narrowPeak 2>> $ol && mv -f *_peaks.narrowPeak.gz $oz 2>> $ol; else gzip < /dev/null > $oz; echo \"File $bam empty\" >> $ol; fi
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
    fw = fn+'_peak.fw.bedg.gz'
    fr = fn+'_peak.re.bedg.gz'
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
    
    PAIRSAMPLES = PEAKSAMPLES.collect{
        element -> return "${workflow.workDir}/../"+element+".bam"
    }

    peaksamples_ch = Channel.fromPath(PAIRSAMPLES.sort())
    genomefile = Channel.fromPath(REF)

    UnzipGenome(genomefile)
    FindPeaks(peaksamples_ch.collate(2))
    PeakToBedg(FindPeaks.out.peak.combine(UnzipGenome.out.chromsize))
    NormalizeBedg(PeakToBedg.out.bedgf.collate(1), PeakToBedg.out.bedgr.collate(1))
    PeakToTRACKS(NormalizeBedg.out.bedgf.collate(1), NormalizeBedg.out.bedgr.collate(1).combine(UnzipGenome.out.chromsize))
    GenerateTrack(PeakToTRACKS.out.bwf.collect(), PeakToTRACKS.out.bwr.collect())

    emit:
    trackdb = GenerateTrack.out.trackdb
}
